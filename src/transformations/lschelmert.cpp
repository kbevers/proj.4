/*****************************************************************************
* Project:	PROJ
* Purpose:	Helmert Least Squared Collocation
* Author:	Sveinung Himle <sveinung.himle at statkart.no>
*
******************************************************************************
* Copyright (c) 2020, Sveinung Himle <sveinung.himle at statkart.no>
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
#include "proj\internal\nlohmann\json.hpp"
#include "Eigen\Eigen"

//using namespace GeoJsonMultiPolygon;
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
		ListOfMultiPolygons polygons{};
	};
}

namespace ns
{
	// TODO: Move to common class
	struct GeoJProperties
	{
		std::string name;
	};

	struct GeoJType
	{
		std::string name;
		GeoJProperties properties;
	};

	struct GeoJCrs
	{
		GeoJType type;
	};

	struct GeoJson
	{
		std::string type;
		std::string name;
		GeoJCrs crs;
		//std::string features;
	};
}
 
struct COMMONPOINTS* find_CommonPointList(projCtx ctx, PJ_LP input, int cp_count, pj_cp **cps)
{
	int iCp;

	for (iCp = 0; iCp < cp_count; iCp++)
	{
		pj_cp *gi = cps[iCp];
		
		COMMONPOINTS *cp = gi->cp;
		/*	if (ct->ll.phi - epsilon > input.phi
			|| ct->ll.lam - epsilon > input.lam
			|| (ct->ll.phi + (ct->lim.phi - 1) * ct->del.phi + epsilon < input.phi)
			|| (ct->ll.lam + (ct->lim.lam - 1) * ct->del.lam + epsilon < input.lam)) {
			continue;
		}*/

		while (gi->child)
		{
			pj_cp *child;

			for (child = gi->child; child != nullptr; child = child->next)
			{
				COMMONPOINTS *cp1 = child->cp;
				/*
				epsilon = (fabs(ct1->del.phi)+fabs(ct1->del.lam))/10000.0;

				if( ct1->ll.phi - epsilon > input.phi
					|| ct1->ll.lam - epsilon > input.lam
					|| (ct1->ll.phi+(ct1->lim.phi-1)*ct1->del.phi + epsilon < input.phi)
					|| (ct1->ll.lam+(ct1->lim.lam-1)*ct1->del.lam + epsilon < input.lam) ) {
					continue;
				}*/
				break;
				
			} 
			if (child == nullptr) 
				break;

			gi = child;
			cp = child->cp;
		}
		if (cp->pJ_LPZ_PairList == nullptr || cp->noOfPoints == 0)
		{
			if (!pj_cp_load(ctx, gi))
			{
				pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_CPL);
				return nullptr;
			}
		}
	}
	return nullptr;
};

MatrixXd CovarianceNN(PJ_LP *lp, std::vector<PJ_LPZ_Pair> *commonPointList, PJ_DIRECTION direction, double k = 0.00039, double c = 0.06900 * M_PI / 180.0)
{
	auto np = commonPointList->size();

	double coslat = cos(lp->phi);

	MatrixXd cnn(np, np);
	
	for (int i = 0; i < np; i++)
	{
		PJ_LPZ_Pair pp1 = commonPointList->at(i);

		PJ_LPZ p1 = (direction == PJ_FWD) ? pp1.fromPoint : pp1.toPoint;

		for (int j = 0; j < np; j++)
		{
			PJ_LPZ_Pair pp2 = commonPointList->at(j);
			PJ_LPZ p2 = (direction == PJ_FWD) ? pp2.fromPoint : pp2.toPoint;

			double dist = hypot(p1.phi - p2.phi, (p1.lam * coslat) - (p2.lam * coslat));
			double a = (M_PI / 2.0) * (dist / c);
			cnn(i, j) = k * exp(-a) * cos(a);
		}
	}
	return cnn;
}

MatrixXd CovarianceMN(PJ_LP *lp, std::vector<PJ_LPZ_Pair> *commonPointList, PJ_DIRECTION direction, double k = 0.00039, double c = 0.06900 * M_PI / 180.0)
{
	auto np = commonPointList->size();
	
	double coslat = cos(lp->phi);

	double x = lp->phi;
	double y = lp->lam * coslat;

	MatrixXd cmn(np, 1);
	 
	for (int i = 0; i < np; i++)
	{
		PJ_LPZ_Pair pointPair = commonPointList->at(i);
		PJ_LPZ pointFrom = (direction == PJ_FWD) ? pointPair.fromPoint : pointPair.toPoint;

		double dist = hypot(pointFrom.phi - x, (pointFrom.lam * coslat) - y);
		double a = (M_PI / 2.0) * (dist / c);
		cmn(i, 0) = k * exp(-a) * cos(a);
	}
	return cmn;
}

/******************************************************************************************
* http://www.mygeodesy.id.au/documents/Coord%20Transforms%20in%20Cadastral%20Surveying.pdf
* https://www.degruyter.com/downloadpdf/j/rgg.2014.97.issue-1/rgg-2014-0009/rgg-2014-0009.pdf
/******************************************************************************************/

static PJ* calculateHelmertParameter(PJ *P, PJ_LP *lp, std::vector<PJ_LPZ_Pair> *commonPointList, PJ_DIRECTION direction)
{
	struct pj_opaque_lschelmert *Q = (struct pj_opaque_lschelmert *) P->opaque;

	auto np = commonPointList->size();

	if (np < 3)
	{
		proj_log_error(P, "lschelemert: common points are less than 3.");		 
		return nullptr;
	}
	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
	{
		proj_log_trace(P, "Input phi, lam: (%12.10f, %12.10f)", lp->phi, lp->lam);
	}

	double coslat = cos(lp->phi);

	// TODO: Include in proj string
	// Covariance matrices:
	double k = 0.00039;
	double c = 0.3;

	MatrixXd cnn = CovarianceNN(lp, commonPointList, direction, k, c);
	MatrixXd cmn = CovarianceMN(lp, commonPointList, direction, k, c);	 
 
	// Vector From System:
	MatrixXd xF(np, 1); MatrixXd yF(np, 1);

	// Vector To System:
	MatrixXd xT(np, 1); MatrixXd yT(np, 1);

	for (int i = 0; i < np; i++)
	{
		PJ_LPZ_Pair pointPair = commonPointList->at(i);

		PJ_LPZ pointFrom = (direction == PJ_FWD) ? pointPair.fromPoint : pointPair.toPoint;
		PJ_LPZ pointTo = (direction == PJ_FWD) ? pointPair.toPoint : pointPair.fromPoint;

		xF(i, 0) = pointFrom.phi;
		yF(i, 0) = pointFrom.lam * coslat;

		xT(i, 0) = pointTo.phi;
		yT(i, 0) = pointTo.lam * coslat;

		if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
			proj_log_trace(P, "Source: (%10.8f, %10.8f) Target: (%10.8f, %10.8f)", pointFrom.phi, pointFrom.lam, pointTo.phi, pointTo.lam);			
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
	
	// Signal of commom points
	MatrixXd snx(np, 1); MatrixXd sny(np, 1);

	for (int i = 0; i < np; i++)
	{
		snx(i) = dxT(i) - a * dxF(i) - b * dyF(i);
		sny(i) = dyT(i) + b * dxF(i) - a * dyF(i);
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

	return P;
} 

bool DistanceLess(const PJ_LPZ_Pair& lhs, const PJ_LPZ_Pair& rhs)
{ 
	return lhs.dist < rhs.dist;
}

/***********************************************************************
* https://stackoverflow.com/questions/4509798/finding-nearest-point-in-an-efficient-way
/***********************************************************************/
std::vector<PJ_LPZ_Pair> findClosestPoints(COMMONPOINTS *commonPointList, PJ_LP lp, int areaId, PJ_DIRECTION direction, int n = 20)
{
	std::vector<PJ_LPZ_Pair> distances;
	std::vector<PJ_LPZ_Pair> closestDistances;

	auto np = commonPointList->pJ_LPZ_PairList->size();

	double coslat = cos(lp.phi);

	for (int i = 0; i < np; i++)
	{
		PJ_LPZ_Pair pair = commonPointList->pJ_LPZ_PairList->at(i);
		PJ_LPZ point = (direction == PJ_FWD) ? pair.fromPoint : pair.toPoint;

		double deltaPhi = point.phi - lp.phi;
		double deltaLam = (point.lam - lp.lam) * coslat;

		pair.dist = hypot(deltaPhi, deltaLam);
		 
		distances.push_back(pair);
	}

	std::sort(distances.begin(), distances.end(), DistanceLess);

	for (int i = 0; i < n; i++)
		closestDistances.push_back(distances[i]);

	return closestDistances;
}

int proj_cp_init(PJ* P, const char *cps)
{
	char *scps = (char *)pj_malloc((strlen(cps) + 1 + 1) * sizeof(char));
	sprintf(scps, "%s%s", "s", cps);

	if (P->cplist == nullptr)
	{
		P->cplist = pj_cplist(P->ctx, pj_param(P->ctx, P->params, scps).s, &(P->cplist_count));

		if (P->cplist == nullptr || P->cplist_count == 0)
		{
			pj_dealloc(scps);
			return 0;
		}
	}

	if (P->cplist_count == 0)
		proj_errno_set(P, PJD_ERR_FAILED_TO_LOAD_CPL);

	pj_dealloc(scps);
	return P->cplist_count;
}

struct COMMONPOINTS* find_cp(projCtx ctx, PJ_LP input, int cp_count, PJ_COMMONPOINTS **cps)
{
	int icp;

	for (icp = 0; icp < cp_count; icp++)
	{
		PJ_COMMONPOINTS *gi = cps[icp];
		struct COMMONPOINTS *cp = gi->cp;	
		
		while (gi->child)
		{
			PJ_COMMONPOINTS *child;

			for (child = gi->child; child != nullptr; child = child->next)
			{
				struct COMMONPOINTS *cp1 = child->cp;
				break;
			}
			if (child == nullptr)
				break;

			gi = child;
			cp = child->cp;
		}
		if (cp->pJ_LPZ_PairList == nullptr)
		{
			if (!pj_cp_load(ctx, gi))
			{
				pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_GRID);
				return nullptr;
			}
		}
		return cp;
	}
	return nullptr;
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

	struct COMMONPOINTS *cp;
	cp = find_cp(P->ctx, point.lp, P->cplist_count, P->cplist);

	if (cp == nullptr)
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GRID);
		return point.xyz;
	}

	int areaId = 1; // areaIdPoint(&point.lp);
	double n = 8;
	auto closestPoints = findClosestPoints(cp, point.lp, areaId, PJ_FWD, n);
	
	if (closestPoints.size() == 0)
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GRID);
		return point.xyz;
	}
	if (!calculateHelmertParameter(P, &point.lp, &closestPoints, PJ_FWD))
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GRID);
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

	struct COMMONPOINTS *cp;
	cp = find_cp(P->ctx, point.lp, P->cplist_count, P->cplist);

	if (cp == nullptr)
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GRID);
		return point.lpz;
	} 

	int areaId = 1; //areaIdPoint(&point.lp);
	double n = 8; // TODO: Make as parameter...
	auto closestPoints = findClosestPoints(cp, point.lp, areaId, PJ_INV, n);

	if (closestPoints.size() == 0)
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GRID);
		return point.lpz;
	}
	if (!calculateHelmertParameter(P, &point.lp, &closestPoints, PJ_INV))
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GRID);
		return point.lpz;
	}
	point.lp = helmert_apply(P, point.lp);
	point.lp = collocation_apply(P, point.lp);

	return point.lpz;
}

static PJ *destructor(PJ *P, int errlev) {
	if (nullptr == P)
		return nullptr;

	auto Q = static_cast<struct pj_opaque_lschelmert*>(P->opaque);
	if (Q)
	{
		//if (Q->->cart)
		//	Q->cart->destructor(Q->cart, errlev);

		delete Q;
	}
	P->opaque = nullptr;

	return pj_default_destructor(P, errlev);
}

PJ *TRANSFORMATION(lschelmert, 0)
{	
	//struct pj_opaque_lschelmert *Q = static_cast<struct pj_opaque_lschelmert*>(pj_calloc(1, sizeof(struct pj_opaque_lschelmert)));
	
	auto Q = new pj_opaque_lschelmert;
	P->opaque = (void *)Q;
	P->destructor = destructor;

	if (Q == nullptr )
	{
		return pj_default_destructor(P, ENOMEM);
	}

	int has_polygons = pj_param(P->ctx, P->params, "tpolygons").i;
	if (has_polygons == 0)
	{
		proj_log_error(P, "cp_trans: +polygon parameter missing.");
		return pj_default_destructor(P, PJD_ERR_NO_ARGS);
	}

	P->opaque = (void *)Q;
	P->fwd4d = nullptr;
	P->inv4d = nullptr;
	P->fwd3d = forward_3d;
	P->inv3d = reverse_3d;
	P->fwd = nullptr;
	P->inv = nullptr;

	// ?
	P->left = PJ_IO_UNITS_RADIANS; 
	P->right = PJ_IO_UNITS_RADIANS;

	if (0 == pj_param(P->ctx, P->params, "tcp_trans").i) 
	{
		proj_log_error(P, "cp_trans: +cp_trans parameter missing.");
		return pj_default_destructor(P, PJD_ERR_NO_ARGS);
	}
	
	proj_cp_init(P, "cp_trans");

	if (proj_errno(P))
	{
		proj_log_error(P, "cp_trans: could not find required cp_tran(s).");
		return pj_default_destructor(P, PJD_ERR_FAILED_TO_LOAD_CPL);
	}
	pj_polygon_init(P, "polygons");
 
	return P;
}
