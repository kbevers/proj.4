/*****************************************************************************

Documentation

******************************************************************************

******************************************************************************
* Project:	PROJ
* Purpose:	Helmert Least Squared Collocation
* Author:	Sveinung Himle <sveinung.himle at kartverket.no>
*
******************************************************************************
* Copyright (c) 2020, Sveinung Himle <sveinung.himle at kartverket.no>
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
******************************************************************************/
#define PJ_LIB__

#include <errno.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstddef>
#include <algorithm>

#include "proj.h"
#include "proj_internal.h"
#include "point_in_polygon.h"
#include "geojsonPolygon.hpp"
#include "pps.hpp"
#include "pps_set.hpp"
#include "proj\internal\nlohmann\json.hpp"
#include <eigen3/Eigen/Dense>

using namespace NS_PROJ;
using namespace Eigen;
using json = nlohmann::json;

PROJ_HEAD(lschelmert, "2D Helmert parameter estimation with collocation");

namespace
{
	struct pj_opaque_lschelmert
	{
		/************************************************************************
			Projection specific elements for the "lschelmert" PJ object
		************************************************************************/
		PJ_XYZ xyz;
		PJ_XYZ xyz_0;
		PJ_XYZ dxyz;
		PJ_XYZ refp;
		PJ_OPK opk;
		PJ_OPK opk_0;
		PJ_OPK dopk;
		double a;
		double b;
		double tx;
		double ty;
		double xA0;
		double yA0;
		double xB0;
		double yB0;
		double signalx;
		double signaly;
		double sigmaHelmert;
		int n_points;
		double maximum_dist;
		double ccoll;
		double kcoll;

		ListOfMultiPolygons polygons{};
		ListOfPpSet pps{};
	};
}

MatrixXd CovarianceNN(PJ_LP *lp, const std::vector<LPZ_Pair> *pairList, PJ_DIRECTION direction, double k = 0.00039, double c = 7.7)
{ 
	int i = 0;
	int j = 0;
	auto np = pairList->size();
	double coslat = cos(lp->phi);
	
	// Scaled to radians
	c *= 2.0 / (M_PI * 6390.0);

	MatrixXd cnn(np, np);

	for (auto&& pair1 : *pairList)
	{  
		j = 0;

		PJ_LPZ p1 = (direction == PJ_FWD) ? pair1.FromPoint() : pair1.ToPoint();

		for (auto&& pair2 : *pairList)
		{ 
			PJ_LPZ p2 = (direction == PJ_FWD) ? pair2.FromPoint() : pair2.ToPoint();

			double dist = hypot(p1.phi - p2.phi, (p1.lam * coslat) - (p2.lam * coslat));
			double a = (dist / c);

			cnn(i, j++) = k * exp(-a) * cos(a);
		}
		i++;
	}
	return cnn;
}

MatrixXd CovarianceMN(PJ_LP *lp, std::vector<LPZ_Pair> *pairList, PJ_DIRECTION direction, double k = 0.00039, double c = 7.7)
{	
	int i = 0;	
	double coslat = cos(lp->phi);
	double x = lp->phi;
	double y = lp->lam * coslat;
	auto np = pairList->size();
	
	// Scaled to radians
	c *= 2.0 / (M_PI * 6390.0);

	MatrixXd cmn(np, 1);

	for (auto&& pair : *pairList)
	{   
		PJ_LPZ p = (direction == PJ_FWD) ? pair.FromPoint() : pair.ToPoint();
		
		double dist = hypot(p.phi - x, (p.lam * coslat) - y);
		double a = (dist / c);

		cmn(i++, 0) = k * exp(-a) * cos(a);
	}
	return cmn;
}

/******************************************************************************************
* http://www.mygeodesy.id.au/documents/Coord%20Transforms%20in%20Cadastral%20Surveying.pdf
* https://www.degruyter.com/downloadpdf/j/rgg.2014.97.issue-1/rgg-2014-0009/rgg-2014-0009.pdf
*https://elib.uni-stuttgart.de/bitstream/11682/10584/1/MscThesis_Dalu_Dong_GIS_04_09_2019.pdf
/******************************************************************************************/
static PJ* calculateHelmertParameter(PJ *P, PJ_LP *lp, std::vector<LPZ_Pair> *pairList, PJ_DIRECTION direction)
{ 
	struct pj_opaque_lschelmert *Q = (struct pj_opaque_lschelmert *) P->opaque;

	auto np = pairList->size();

	if (np < 3)
	{
		proj_log_error(P, "lschelmert: common point pairs are less than 3.");		 
		return nullptr;
	}
	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
		proj_log_trace(P, "Input phi, lam: (%12.10f, %12.10f)", lp->phi, lp->lam);	 

	double coslat = cos(lp->phi);
	double k = Q->kcoll == HUGE_VAL ? 0.00039 : Q->kcoll;
	double c = Q->ccoll == HUGE_VAL ? 7.7 : Q->ccoll;

    // Covariance matrices:
	MatrixXd cnn = CovarianceNN(lp, pairList, direction, k, c);
	MatrixXd cmn = CovarianceMN(lp, pairList, direction, k, c);
 
	// Vector source system:
	MatrixXd xA(np, 1); MatrixXd yA(np, 1);
	
	// Vector target system:
	MatrixXd xB(np, 1); MatrixXd yB(np, 1);

	int l = 0;

	for (auto&& pair : *pairList)
	{
		PJ_LPZ pointFrom = (direction == PJ_FWD) ? pair.FromPoint() : pair.ToPoint();
	    PJ_LPZ pointTo = (direction == PJ_FWD) ? pair.ToPoint() : pair.FromPoint();

		xA(l, 0) = pointFrom.phi;
		yA(l, 0) = pointFrom.lam * coslat;

		xB(l, 0) = pointTo.phi;
		yB(l, 0) = pointTo.lam * coslat;

		if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
			proj_log_trace(P, "Source: (%10.8f, %10.8f) Target: (%10.8f, %10.8f)", pointFrom.phi, pointFrom.lam, pointTo.phi, pointTo.lam);
	
		l++;
	}
	
	// Weight matrix p is the inverted cnn
	MatrixXd p = cnn.inverse();

	// W, sum of weight p:	
	auto w_sum = p.sum();

	// Weight for each point:
	ArrayXXd w_points = p * MatrixXd::Ones(np, 1);

	// Transposed w_points
	MatrixXd w_pointsTrans = w_points.transpose();
	 
 	// Mass center:
	double xA0 = (w_pointsTrans * xA / w_sum).value();
	double yA0 = (w_pointsTrans * yA / w_sum).value();
	double xB0 = (w_pointsTrans * xB / w_sum).value();
	double yB0 = (w_pointsTrans * yB / w_sum).value();

	// Coordinates in Mass center origin:
	ArrayXXd xA_ = xA - MatrixXd::Ones(np, 1) * xA0;
	ArrayXXd yA_ = yA - MatrixXd::Ones(np, 1) * yA0;
	ArrayXXd xB_ = xB - MatrixXd::Ones(np, 1) * xB0;
	ArrayXXd yB_ = yB - MatrixXd::Ones(np, 1) * yB0;

	// Normal equation parameters:
	double n = 0.0;
	double t1 = 0.0;
	double t2 = 0.0;

	// Putting values into normal equation
	for (int i = 0; i < np; i++)
	{
		n += (pow(xA_(i), 2) * w_points(i)) + (pow(yA_(i), 2) * w_points(i));
		t1 += (xA_(i) * xB_(i) + yA_(i) * yB_(i)) * w_points(i);
		t2 += (yA_(i) * xB_(i) - xA_(i) * yB_(i)) * w_points(i);
	}

	// Estimated Helmert parameters
	double a = t1 / n;
	double b = t2 / n;
	double tx = xB0 - a * xA0 - b * yA0;
	double ty = yB0 + b * xA0 - a * yA0;

	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
		proj_log_trace(P, "Estimated Helmert parameters a, b, Tx, Ty: (%12.10f, %12.10f, %12.10f, %12.10f)", a, b, tx, ty);
	
	// Signal (residuals) of common points
	MatrixXd snx(np, 1); MatrixXd sny(np, 1);

	// Residuals referred in mass center origin
	for (int i = 0; i < np; i++)
	{
		snx(i) = xB_(i) - a * xA_(i) - b * yA_(i);
		sny(i) = yB_(i) + b * xA_(i) - a * yA_(i);
	}

	double sigma = 0.0;
		
	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
	{
		for (int i = 0; i < np; i++)
		{
			for (int j = 0; j < np; j++)
			{
				sigma += snx(j) * snx(j) * p(j, i);
				sigma += sny(j) * sny(j) * p(j, i);
			}
		}
		sigma /= 2 * np;
		sigma = sqrt(sigma);

		proj_log_trace(P, "Estimated sigma Helmert transformation: (%12.10f radians)", sigma);
	}
		
	// Signal of target point
	double smx = (cmn.transpose() * p * snx).value();
	double smy = (cmn.transpose() * p * sny).value();
	 
	Q->a = a;
	Q->b = b;
	Q->tx = tx;
	Q->ty = ty;
	Q->xA0 = xA0;
	Q->yA0 = yA0;
	Q->xB0 = xB0;
	Q->yB0 = yB0;
	Q->signalx = smx;
	Q->signaly = smy;
	Q->sigmaHelmert = sigma;

	return P;
} 

bool DistanceLess(const LPZ_Pair& lhs, const LPZ_Pair& rhs)
{
	return lhs.Distance() < rhs.Distance();
}

/***********************************************************************
* https://stackoverflow.com/questions/4509798/finding-nearest-point-in-an-efficient-way
/***********************************************************************/
std::vector<LPZ_Pair> findClosestPoints(PointPairs *ppList, PJ_LP lp, __int32 areaId, PJ_DIRECTION direction, int n = 20, double maximum_dist = 0.1)
{
	std::vector<LPZ_Pair> distances {};

	double coslat = cos(lp.phi);

	for (auto pair : ppList->LpzPairList())
	{
		if (areaId != 0 && pair.Area() != areaId)
			continue;

		PJ_LPZ point = (direction == PJ_FWD) ? pair.FromPoint() : pair.ToPoint();

		double deltaPhi = point.phi - lp.phi;
		double deltaLam = (point.lam - lp.lam) * coslat;

		if (hypot(deltaPhi, deltaLam) > maximum_dist)
			continue;

		pair.SetDistance(hypot(deltaPhi, deltaLam));
		
		distances.push_back(pair);
	}	 
	std::sort(distances.begin(), distances.end(), DistanceLess);

	auto np = distances.size();

	if (n < np)
		distances.resize(n);

	return distances;
}

PJ_LP helmert_apply(PJ *P, PJ_LP lp)
{
	struct pj_opaque_lschelmert *Q = (struct pj_opaque_lschelmert *) P->opaque;

	PJ_LP out;
	out.lam = HUGE_VAL; out.phi = HUGE_VAL;

	double coslat = cos(lp.phi);

	double xTrans = Q->xB0 - Q->a * (Q->xA0 - lp.phi) - Q->b * (Q->yA0 - lp.lam * coslat);
	double yTrans = Q->yB0 + Q->b * (Q->xA0 - lp.phi) - Q->a * (Q->yA0 - lp.lam * coslat);
	
	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
		proj_log_trace(P, "Helmert transformated phi, lam: (%12.10f, %12.10f)", xTrans, yTrans / coslat);

	out.phi = xTrans;
	out.lam = yTrans / coslat;

	return out;
}

PJ_LP collocation_apply(PJ *P, PJ_LP lp)
{
	struct pj_opaque_lschelmert *Q = (struct pj_opaque_lschelmert *) P->opaque;

	PJ_LP out;
	out.lam = HUGE_VAL; out.phi = HUGE_VAL;

	double coslat = cos(lp.phi);

	out.phi = lp.phi + Q->signalx;
	out.lam = lp.lam + (Q->signaly / coslat);

	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
		proj_log_trace(P, "LSC predicted phi, lam: (%12.10f, %12.10f)", out.phi, out.lam);

	return out;
}

static PJ_XYZ forward_3d(PJ_LPZ lpz, PJ *P)
{
	struct pj_opaque_lschelmert *Q = (struct pj_opaque_lschelmert *) P->opaque;

	PJ_COORD point = { {0,0,0,0} };
	point.lpz = lpz;

	Q->maximum_dist = Q->maximum_dist == HUGE_VAL ? 0.1 : Q->maximum_dist;

	PointPairs* pointPairs = findPointPairs(Q->pps, lpz, Q->maximum_dist);
	
	if (pointPairs == nullptr)
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return point.xyz;
	}

 	__int32 areaId = areaIdPoint(Q->polygons, &point.lp);
	int n = Q->n_points == FP_NORMAL ? 20 : Q->n_points; // Default 20 point candidates	 
	auto closestPoints = findClosestPoints(pointPairs, point.lp, areaId, PJ_FWD, n, Q->maximum_dist);
	
	if (closestPoints.size() == 0)
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return point.xyz;
	}
	if (!calculateHelmertParameter(P, &point.lp, &closestPoints, PJ_FWD))
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return point.xyz;
	}
	point.lp = helmert_apply(P, point.lp);
	point.lp = collocation_apply(P, point.lp);
	 
	return point.xyz;
}

static PJ_LPZ reverse_3d(PJ_XYZ xyz, PJ *P)
{
	struct pj_opaque_lschelmert *Q = (struct pj_opaque_lschelmert *) P->opaque;

	PJ_COORD point = { {0,0,0,0} };
	point.xyz = xyz;
	auto lpz = point.lpz;

	Q->maximum_dist = Q->maximum_dist == HUGE_VAL ? 0.1 : Q->maximum_dist; // Default 0.1 radians

	PointPairs *pointPairs = findPointPairs(Q->pps, lpz, Q->maximum_dist);

	if (pointPairs == nullptr)
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return point.lpz;
	}

	__int32 areaId = areaIdPoint(Q->polygons, &point.lp);
	int n = Q->n_points == FP_NORMAL ? 20 : Q->n_points; // Default 20 point candidates	
	auto closestPoints = findClosestPoints(pointPairs, point.lp, areaId, PJ_INV, n, Q->maximum_dist);

	if (closestPoints.size() == 0)
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return point.lpz;
	}
	if (!calculateHelmertParameter(P, &point.lp, &closestPoints, PJ_INV))
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return point.lpz;
	}
	point.lp = helmert_apply(P, point.lp);
	point.lp = collocation_apply(P, point.lp);

	return point.lpz;
}

static PJ *destructor(PJ *P, int errlev)
{
	if (nullptr == P)
		return nullptr;

	auto Q = static_cast<struct pj_opaque_lschelmert*>(P->opaque);
	if (Q)
	{
		delete Q;
	}
	P->opaque = nullptr;

	return pj_default_destructor(P, errlev);
}

// TODO: Static method might be deleted
static void reassign_context(PJ* P, PJ_CONTEXT* ctx)
{
	auto Q = (struct pj_opaque_lschelmert *) P->opaque;
	for (auto& poly : Q->polygons)
	{
		poly->reassign_context(ctx);
	}
}

PJ *TRANSFORMATION(lschelmert, 0)
{
	//struct pj_opaque_lschelmert *Q = static_cast<struct pj_opaque_lschelmert*>(pj_calloc(1, sizeof(struct pj_opaque_lschelmert)));
	auto Q = new pj_opaque_lschelmert;
	P->opaque = (void *)Q;
	P->destructor = destructor;
	P->reassign_context = reassign_context;	

	if (Q == nullptr)	 
		return destructor(P, ENOMEM);	

	P->opaque = (void *)Q;
	P->fwd4d = nullptr;
	P->inv4d = nullptr;
	P->fwd3d = forward_3d;
	P->inv3d = reverse_3d;
	P->fwd = nullptr;
	P->inv = nullptr;

	P->left = PJ_IO_UNITS_RADIANS; 
	P->right = PJ_IO_UNITS_RADIANS;

	if (0 == pj_param(P->ctx, P->params, "tpp_trans").i) 
	{
		proj_log_error(P, "pp_trans: +pp_trans parameter missing.");
		return pj_default_destructor(P, PJD_ERR_NO_ARGS);
	}
	Q->pps = pj_cp_init(P, "pp_trans");

	Q->n_points = FP_NORMAL;
	if (pj_param_exists(P->params, "n_points"))	
		Q->n_points = pj_param(P->ctx, P->params, "in_points").i;

	Q->maximum_dist = HUGE_VAL;
	if (pj_param_exists(P->params, "max_dist"))
		Q->maximum_dist = pj_param(P->ctx, P->params, "dmax_dist").f;

	Q->ccoll = HUGE_VAL;
	if (pj_param_exists(P->params, "ccoll"))
		Q->ccoll = pj_param(P->ctx, P->params, "dccoll").f;

	Q->kcoll = HUGE_VAL;
	if (pj_param_exists(P->params, "kcoll"))
		Q->kcoll = pj_param(P->ctx, P->params, "dkcoll").f;

	int has_polygons = pj_param(P->ctx, P->params, "tpolygons").i;
	if (has_polygons > 1)	 
		Q->polygons = pj_polygon_init(P, "polygons");	  

	if (proj_errno(P))
	{		 
		proj_log_error(P, "pair_trans: could not find required pair_trans file.");
		return destructor(P, PJD_ERR_FAILED_TO_LOAD_CPT);
	}
	return P;
}
