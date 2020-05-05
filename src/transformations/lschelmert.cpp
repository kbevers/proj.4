/*****************************************************************************
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
#include "geocent.h"
#include "point_in_polygon.h"
#include "geojsonPolygon.hpp"
#include "pps.hpp"
#include "pps_set.hpp"
#include "proj\internal\nlohmann\json.hpp"
#include "Eigen\Eigen"

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
		double xF0;
		double yF0;
		double xT0;
		double yT0;
		double signalx;
		double signaly;
		double sigmaHelmert;
		int n_points;
		double c_coll;
		double k_coll;

		ListOfMultiPolygons polygons{};
		ListOfPpSet pps{};
	};
}

MatrixXd CovarianceNN(PJ_LP *lp, const std::vector<LPZ_Pair> *pairList, PJ_DIRECTION direction, double k = 0.00039, double c = 0.001204)
{ 
	int i = 0;
	int j = 0;
	auto np = pairList->size();
	double coslat = cos(lp->phi);
	
	MatrixXd cnn(np, np);
	
	for (auto&& pair1 : *pairList)
	{  
		j = 0;

		PJ_LPZ p1 = (direction == PJ_FWD) ? pair1.FromPoint() : pair1.ToPoint();

		for (auto&& pair2 : *pairList)
		{ 
			PJ_LPZ p2 = (direction == PJ_FWD) ? pair2.FromPoint() : pair2.ToPoint();

			double dist = hypot(p1.phi - p2.phi, (p1.lam * coslat) - (p2.lam * coslat));
			double a = (M_PI / 2.0) * (dist / c);
			cnn(i, j++) = k * exp(-a) * cos(a);		
		}
		i++;
	}
	return cnn;
}

MatrixXd CovarianceMN(PJ_LP *lp, std::vector<LPZ_Pair> *pairList, PJ_DIRECTION direction, double k = 0.00039, double c = 0.001204)
{	
	int i = 0;	
	double coslat = cos(lp->phi);
	double x = lp->phi;
	double y = lp->lam * coslat;
	auto np = pairList->size();

	MatrixXd cmn(np, 1);

	for (auto&& pair : *pairList)
	{   
		PJ_LPZ p = (direction == PJ_FWD) ? pair.FromPoint() : pair.ToPoint();
		
		double dist = hypot(p.phi - x, (p.lam * coslat) - y);
		double a = (M_PI / 2.0) * (dist / c);
		cmn(i++, 0) = k * exp(-a) * cos(a); 
	}
	return cmn;
}

/******************************************************************************************
* http://www.mygeodesy.id.au/documents/Coord%20Transforms%20in%20Cadastral%20Surveying.pdf
* https://www.degruyter.com/downloadpdf/j/rgg.2014.97.issue-1/rgg-2014-0009/rgg-2014-0009.pdf
/******************************************************************************************/
static PJ* calculateHelmertParameter(PJ *P, PJ_LP *lp, std::vector<LPZ_Pair> *pairList, PJ_DIRECTION direction)
{ 
	struct pj_opaque_lschelmert *Q = (struct pj_opaque_lschelmert *) P->opaque;

	auto np = pairList->size();

	if (np < 3)
	{
		proj_log_error(P, "lschelemert: common points are less than 3.");		 
		return nullptr;
	}
	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)	 
		proj_log_trace(P, "Input phi, lam: (%12.10f, %12.10f)", lp->phi, lp->lam);	 

	double coslat = cos(lp->phi);
	double k = Q->k_coll == HUGE_VAL ? 0.00039 : Q->k_coll;
	double c = Q->c_coll == HUGE_VAL ? 0.001204 : Q->c_coll;

    // Covariance matrices:
	MatrixXd cnn = CovarianceNN(lp, pairList, direction, k, c);
	MatrixXd cmn = CovarianceMN(lp, pairList, direction, k, c);
 
	// Vector From System:
	MatrixXd xF(np, 1); MatrixXd yF(np, 1);

	// Vector To System:
	MatrixXd xT(np, 1); MatrixXd yT(np, 1);

	int i = 0;

	for (auto&& pair : *pairList)
	{
		PJ_LPZ pointFrom = (direction == PJ_FWD) ? pair.FromPoint() : pair.ToPoint();
	    PJ_LPZ pointTo = (direction == PJ_FWD) ? pair.ToPoint() : pair.FromPoint();

		xF(i, 0) = pointFrom.phi;
		yF(i, 0) = pointFrom.lam * coslat;

		xT(i, 0) = pointTo.phi;
		yT(i, 0) = pointTo.lam * coslat;

		if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
			proj_log_trace(P, "Source: (%10.8f, %10.8f) Target: (%10.8f, %10.8f)", pointFrom.phi, pointFrom.lam, pointTo.phi, pointTo.lam);
	
		i++;
	}

	MatrixXd p = cnn.inverse();

	// N, sum of p:	
	auto sump = p.sum();

	// sep, row sum of p:	
	ArrayXXd sep = p * MatrixXd::Ones(np, 1);

	// Mass center:
	double xF0 = (MatrixXd::Ones(1, np) * p * xF / sump).value();
	double yF0 = (MatrixXd::Ones(1, np) * p * yF / sump).value();
	double xT0 = (MatrixXd::Ones(1, np) * p * xT / sump).value();
	double yT0 = (MatrixXd::Ones(1, np) * p * yT / sump).value();

	// Coordinates in Mass center origin:
	ArrayXXd dxF = xF - MatrixXd::Ones(np, 1) * xF0;
	ArrayXXd dyF = yF - MatrixXd::Ones(np, 1) * yF0;
	ArrayXXd dxT = xT - MatrixXd::Ones(np, 1) * xT0;
	ArrayXXd dyT = yT - MatrixXd::Ones(np, 1) * yT0;

	// Normal equation:
	double n = 0.0;
	double t1 = 0.0;
	double t2 = 0.0;

	for (int i = 0; i < np; i++)
	{
		n += (pow(dxF(i), 2) * sep(i)) + (pow(dyF(i), 2) * sep(i));
		t1 += (dxF(i) * dxT(i) + dyF(i) * dyT(i)) * sep(i);
		t2 += (dyF(i) * dxT(i) - dxF(i) * dyT(i)) * sep(i);
	}

	double a = t1 / n;
	double b = t2 / n;
	double tx = xT0 - a * xF0 - b * yF0;
	double ty = yT0 + b * xF0 - a * yF0;

	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
		proj_log_trace(P, "Estimated Helmert parameters a, b, Tx, Ty: (%12.10f, %12.10f, %12.10f, %12.10f)", a, b, tx, ty);
	
	// Signal (residuals) of common points
	MatrixXd snx(np, 1); MatrixXd sny(np, 1);

	for (int i = 0; i < np; i++)
	{
		snx(i) = dxT(i) - a * dxF(i) - b * dyF(i);
		sny(i) = dyT(i) + b * dxF(i) - a * dyF(i);
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
	Q->xF0 = xF0;
	Q->yF0 = yF0;
	Q->xT0 = xT0;
	Q->yT0 = yT0;
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

std::vector<LPZ_Pair> findClosestPoints(PointPairs *ppList, PJ_LP lp, __int32 areaId, PJ_DIRECTION direction, int n = 20)
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

		if (hypot(deltaPhi, deltaLam) > 1.0)
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

	double xTrans = Q->xT0 - Q->a * (Q->xF0 - lp.phi) - Q->b * (Q->yF0 - lp.lam * coslat);
	double yTrans = Q->yT0 + Q->b * (Q->xF0 - lp.phi) - Q->a * (Q->yF0 - lp.lam * coslat);
	
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

	PointPairs* pp = findPointPairs(Q->pps, lpz);

	if (pp == nullptr)
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return point.xyz;
	}

 	__int32 areaId = areaIdPoint(Q->polygons, &point.lp);
	int n = Q->n_points == FP_NORMAL ? 20 : Q->n_points; // Default 20 point candidates
	auto closestPoints = findClosestPoints(pp, point.lp, areaId, PJ_FWD, n);
	
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

	PointPairs *pp = findPointPairs(Q->pps, lpz);

	if (pp == nullptr)
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return point.lpz;
	}

	__int32 areaId = areaIdPoint(Q->polygons, &point.lp);
	int n = Q->n_points == FP_NORMAL ? 20 : Q->n_points;
	auto closestPoints = findClosestPoints(pp, point.lp, areaId, PJ_INV, n);

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

	Q->c_coll = HUGE_VAL;
	if (pj_param_exists(P->params, "c_coll"))
		Q->c_coll = pj_param(P->ctx, P->params, "dc_coll").f;

	Q->k_coll = HUGE_VAL;
	if (pj_param_exists(P->params, "k_coll"))
		Q->k_coll = pj_param(P->ctx, P->params, "dk_coll").f;

	int has_polygons = pj_param(P->ctx, P->params, "tpolygons").i;
	if (has_polygons == 0)
	{
		//proj_log_error(P, "pair_trans: +polygon parameter missing.");
		//return destructor(P, PJD_ERR_NO_ARGS);
	}
	else
		Q->polygons = pj_polygon_init(P, "polygons");

	if (proj_errno(P))
	{
		// TODO: Feil meldingstekst her...
		proj_log_error(P, "pair_trans: could not find required pair_trans file.");
		return destructor(P, PJD_ERR_FAILED_TO_LOAD_CPT);		 
	}
	return P;
}
