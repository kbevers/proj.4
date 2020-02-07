/************************************************************************
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
***********************************************************************/
#define PJ_LIB__

#include <errno.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstddef>
#include <algorithm>

#include "proj_internal.h"
#include "geocent.h"
#include "point_in_polygon.h"
#include "proj\internal\nlohmann\json.hpp"
#include "Eigen\Eigen"

using namespace Eigen;
using json = nlohmann::json;

PROJ_HEAD(lschelmert, "2D Helmert parameter estimation with collocation");

static PJ_XYZ lschelmert_forward_3d(PJ_LPZ lpz, PJ *P);
static PJ_LPZ lschelmert_reverse_3d(PJ_XYZ xyz, PJ *P);

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
		if (cp->pJ_LP_PairList == nullptr || cp->noOfPoints == 0)
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

template<class UnaryFunction>
void recursive_iterate(const json& j, UnaryFunction f)
{
	for (auto it = j.begin(); it != j.end(); ++it)
	{
		// auto typenam = it->type_name();
		auto v = it.value();

		if (it->is_array() || it->is_object())
		{
			recursive_iterate(*it, f);
		}
		else if (it->is_null())
		{
			f(it);
		}
		else if (it->is_number_float())
		{
			auto value = it.value();
		}
		else
		{
			f(it);
		}
	}
}

// TODO: Move to a persistant class
/***********************************************************************
*
/***********************************************************************/
static void testReadGeojson(/*char* fileName*/)
{
	char* fileName = "C:/Prosjekter/SkTrans/EurefNgo/Punksky_tilfeldig/Ngo_areas.geojson";

	FILE *f = fopen(fileName, "rb");
	fseek(f, 0, SEEK_END);
	long fsize = ftell(f);
	fseek(f, 0, SEEK_SET);

	char *string = (char *)malloc(fsize + 1);
	fread(string, fsize, 1, f);
	fclose(f);

	string[fsize] = 0;

	// parse and serialize JSON
	json j_complete = json::parse(string);
	// std::cout << std::setw(4) << j_complete << "\n\n";

	recursive_iterate(j_complete, [](json::const_iterator it)
	{});

	// Test
	auto feat = j_complete.at("features");
	//auto c = j_complete.find("coordinates");

	for (auto it = feat.begin(); it != feat.end(); ++it)
	{
		//auto hh =	it1.key["geometry"];

		auto geo = (*it)["geometry"];
		auto coords = geo.at("coordinates");

		//recursive_iterateToFloat(*it, coords.value);

		//if (it1->is_number_float)
		/*
		for (auto it2 = geo.begin(); it2 != geo.end(); ++it2)
		{
			auto coord = (*it2)["coordinates"];
		}*/
	}

	for (auto& x : feat.items())
	{
		//std::cout << "key: " << x.key() << ", value: " << x.value() << '\n';
	}
	/*
	 json::iterator it = feat.begin();

	 std::for_each(feat.begin(), feat.end(), [](std::string &string)
	 {
		 std::cout << string << "\n";
	 }
	 );
*/
//feat.at()
	auto feat1 = feat.at(0);
	auto geo = feat1.at("geometry");
	auto coords1 = geo.at("coordinates");
	auto coords2 = coords1.at(0);
	auto coords3 = coords2.at(0);
	auto coords4 = coords3.at(0);
	auto coords5 = coords4.at(0);

	//auto multiPolygon = j_complete.array("MultiPolygon");
	//auto items = j_complete.items();
	auto json_string = j_complete.dump();
}

MatrixXd CovarianceNN(PJ_LP *lp, std::vector<PJ_LP_Pair> *commonPointList, PJ_DIRECTION direction, double k = 0.00039, double c = 0.06900 * M_PI / 180.0)
{
	auto np = commonPointList->size();

	double coslat = cos(lp->phi);

	MatrixXd cnn(np, np);
	
	for (int i = 0; i < np; i++)
	{
		PJ_LP_Pair pp1 = commonPointList->at(i);

		PJ_LP p1 = (direction == PJ_FWD) ? pp1.fromPoint : pp1.toPoint;

		for (int j = 0; j < np; j++)
		{
			PJ_LP_Pair pp2 = commonPointList->at(j);
			PJ_LP p2 = (direction == PJ_FWD) ? pp2.fromPoint : pp2.toPoint;

			double dist = hypot(p1.phi - p2.phi, (p1.lam * coslat) - (p2.lam * coslat));
			double a = (M_PI / 2.0) * (dist / c);
			cnn(i, j) = k * exp(-a) * cos(a);
		}
	}
	return cnn;
}

MatrixXd CovarianceMN(PJ_LP *lp, std::vector<PJ_LP_Pair> *commonPointList, PJ_DIRECTION direction, double k = 0.00039, double c = 0.06900 * M_PI / 180.0)
{
	auto np = commonPointList->size();
	
	double coslat = cos(lp->phi);

	double x = lp->phi;
	double y = lp->lam * coslat;

	MatrixXd cmn(np, 1);
	 
	for (int i = 0; i < np; i++)
	{
		PJ_LP_Pair pointPair = commonPointList->at(i);
		PJ_LP pointFrom = (direction == PJ_FWD) ? pointPair.fromPoint : pointPair.toPoint;

		double dist = hypot(pointFrom.phi - x, (pointFrom.lam * coslat) - y);
		double a = (M_PI / 2.0) * (dist / c);
		cmn(i, 0) = k * exp(-a) * cos(a);
	}
	return cmn;
} 

static bool calculateHelmertParameter(PJ *P, PJ_LP *lp, std::vector<PJ_LP_Pair> *commonPointList, PJ_DIRECTION direction)
{	
	struct pj_opaque_lschelmert *Q = (struct pj_opaque_lschelmert *) P->opaque;

	auto np = commonPointList->size();
	bool retvalue = false;

	if (np < 3)
	{
		proj_log_error(P, "lschelemert: common points are less than 3.");
		return retvalue;
	}
	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
	{
		proj_log_trace(P, "Input phi, lam: (%12.10f, %12.10f)", lp->phi, lp->lam);
	}

	double coslat = cos(lp->phi);

	// Covariance matrices:
	MatrixXd cnn = CovarianceNN(lp, commonPointList, direction);
	MatrixXd cmn = CovarianceMN(lp, commonPointList, direction);	 
 
	// Vector From System:
	MatrixXd xF(np, 1); MatrixXd yF(np, 1);

	// Vector To System:
	MatrixXd xT(np, 1); MatrixXd yT(np, 1);

	for (int i = 0; i < np; i++)
	{
		PJ_LP_Pair pointPair = commonPointList->at(i);

		PJ_LP pointFrom = (direction == PJ_FWD) ? pointPair.fromPoint : pointPair.toPoint;
		PJ_LP pointTo = (direction == PJ_FWD) ? pointPair.toPoint : pointPair.fromPoint;

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

	retvalue = true;
	
	return retvalue;
} 
 
/******************************************************************************************
* http://www.mygeodesy.id.au/documents/Coord%20Transforms%20in%20Cadastral%20Surveying.pdf
* https://www.degruyter.com/downloadpdf/j/rgg.2014.97.issue-1/rgg-2014-0009/rgg-2014-0009.pdf
/******************************************************************************************/
static void calculateLscHelmert(PJ *P, PJ_LP *lp, std::vector<PJ_LP_Pair> *commonPointList, PJ_DIRECTION direction)
{
	struct pj_opaque_lschelmert *Q = (struct pj_opaque_lschelmert *) P->opaque;

	auto np = commonPointList->size();

	double k = 0.00039;
	double c = 0.06900 * M_PI / 180.0;

	double coslat = cos(lp->phi);

	double x = lp->phi;
	double y = lp->lam * coslat;

	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
	{
		proj_log_trace(P, "Covariance function coefficients: k = %10.8f, c = %10.8f", k, c);
		proj_log_trace(P, "Input phi, lam: (%12.10f, %12.10f)", lp->phi, lp->lam);
	}

	// Covariance matrices:
	MatrixXd cnn(np, np);
	MatrixXd cmn(np, 1);

	//cout << "cmn:" << endl << cmn << endl;
	//cout << "cnn:" << endl << cnn << endl;

	// Vector From System:
	MatrixXd xF(np, 1); MatrixXd yF(np, 1);

	// Vector To System:
	MatrixXd xT(np, 1); MatrixXd yT(np, 1);
	 
	for (int i = 0; i < np; i++)
	{
		PJ_LP_Pair pointPair = commonPointList->at(i);
				 
		PJ_LP pointFrom = (direction == PJ_FWD) ? pointPair.fromPoint : pointPair.toPoint;
		PJ_LP pointTo = (direction == PJ_FWD) ? pointPair.toPoint : pointPair.fromPoint;

		xF(i, 0) = pointFrom.phi;
		yF(i, 0) = pointFrom.lam * coslat;

		xT(i, 0) = pointTo.phi;
		yT(i, 0) = pointTo.lam * coslat;
		
		if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
		 	proj_log_trace(P, "From phi, lam: (%10.8f, %10.8f) To phi, lam: (%10.8f, %10.8f)", pointFrom.phi, pointFrom.lam, pointTo.phi, pointTo.lam);
	}	 
	for (int i = 0; i < np; i++)
	{
		double dist = hypot(xF(i, 0) - x, yF(i, 0) - y);
		double a = (M_PI / 2.0) * (dist / c);
		cmn(i, 0) = k * exp(-a) * cos(a);

		for (int j = 0; j < np; j++)
		{
			dist = hypot(xF(i, 0) - xF(j, 0), yF(i, 0) - yF(j, 0));
			a = (M_PI / 2.0) * (dist / c);
			cnn(i, j) = k * exp(-a) * cos(a);
		}
	}

	MatrixXd p = cnn.inverse();
	//cout << "p:" << endl << p << endl;

	// N, sum of p:	
	auto sump = p.sum();

	// sep, row sum of p:	
	ArrayXXd sep = p * MatrixXd::Ones(np, 1);
	//cout << sep.square() << endl;

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
	
	//cout << "n:" << endl << n << endl;

	double a = t1 / n;
	double b = t2 / n;

	double tx = xT0 - a * xF0 - b * yF0;
	double ty = yT0 + b * xF0 - a * yF0;

	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)	
		proj_log_trace(P, "Estimated Helmert parameters a, b, Tx, Ty: (%12.10f, %12.10f, %12.10f, %12.10f)", a, b, tx, ty);
	 
	Q->a = a;
	Q->b = b;
	Q->tx = tx;
	Q->ty = ty;	
	Q->xF0 = xF0;
	Q->yF0 = yF0;
	Q->xT0 = xT0;
	Q->yT0 = yT0;

	// Sigma noise
	MatrixXd snx(np, 1); MatrixXd sny(np, 1);

	for (int i = 0; i < np; i++)
	{
		snx(i) = dxT(i) - a * dxF(i) - b * dyF(i);
		sny(i) = dyT(i) + b * dxF(i) - a * dyF(i);
	}
	// cout << "snx:" << endl << snx << endl;
	// cout << "sny:" << endl << sny << endl;

	// (xT0-xT(0, 0))-a*(xF0-xF(0, 0))-b*(yF0-yF(0,0));
	// (yT0-yT(0, 0))+b*(xF0-xF(0, 0))-a*(yF0-yF(0,0));

	double xTrans = xT0 - a * (xF0 - x) - b * (yF0 - y);
	double yTrans = yT0 + b * (xF0 - x) - a * (yF0 - y);

	lp->phi = xTrans;
	lp->lam = yTrans / coslat;

	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
		proj_log_trace(P, "Transformated phi, lam: (%12.10f, %12.10f)", xTrans, yTrans / coslat);

	double smx = (cmn.transpose() * p * snx).value();
	double smy = (cmn.transpose() * p * sny).value();
	// cout << "smx:" << endl << smx << endl;
	// cout << "smy:" << endl << smy << endl;

	lp->phi = xTrans + smx;
	lp->lam = (yTrans + smy) / coslat;

	if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
		proj_log_trace(P, "LSC predicted phi, lam: (%12.10f, %12.10f)", lp->phi, lp->lam);
}

bool DistanceLess(const PJ_LP_Pair& lhs, const PJ_LP_Pair& rhs)
{ 
	return lhs.dist < rhs.dist;
}

/***********************************************************************
* https://stackoverflow.com/questions/4509798/finding-nearest-point-in-an-efficient-way
/***********************************************************************/
std::vector<PJ_LP_Pair> findClosestPoints(COMMONPOINTS *commonPointList, PJ_LP lp, int areaId, PJ_DIRECTION direction, int n = 20)
{
	std::vector<PJ_LP_Pair> distances;
	std::vector<PJ_LP_Pair> closestDistances;

	auto np = commonPointList->pJ_LP_PairList->size();

	double coslat = cos(lp.phi);

	for (int i = 0; i < np; i++)
	{
		PJ_LP_Pair pair = commonPointList->pJ_LP_PairList->at(i);
		PJ_LP point = (direction == PJ_FWD) ? pair.fromPoint : pair.toPoint;

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

bool pointIsInArea(PJ_LP pointPJ_LP, char* fileName)
{
	std::ifstream file(fileName, std::ios::in);

	Point points[] = { {} };
	Point point = { pointPJ_LP.phi, pointPJ_LP.lam };
	vector<Point> pointVector;

	if (file.is_open())
	{
		std::vector<std::vector<std::string>> dataList;
		std::string line = "";
		double x, y;

		while (!file.eof())
		{
			getline(file, line, ',');
			x = atof(line.c_str());

			getline(file, line, '\n');
			y = atof(line.c_str());

			Point areaPoint;
			areaPoint.x = x;
			areaPoint.y = y;

			pointVector.push_back(areaPoint);
		}
		file.close();
	};

	Point *vectorPointer = pointVector.data();
	int n = size(pointVector);

	return isInside(vectorPointer, n, point);
}

int areaIdPoint(PJ_LP *lp) // TODO: Endre namn og argument 
{
	// TODO: Flytte områdefilene
	char* fileName2 = "C:/Prosjekter/SkTrans/EurefNgo/Punksky_tilfeldig/Area2.csv";
	char* fileName3 = "C:/Prosjekter/SkTrans/EurefNgo/Punksky_tilfeldig/Area3.csv";
	char* fileName4 = "C:/Prosjekter/SkTrans/EurefNgo/Punksky_tilfeldig/Area4.csv";

	// TODO: Områdefil som geojson. Json ligg under include/proj/internal/nlohmann
	//testReadGeojson();

	if (pointIsInArea(*lp, fileName2))
		return 2;
	else if (pointIsInArea(*lp, fileName3))
		return 3;
	else if (pointIsInArea(*lp, fileName4))
		return 4;

	return 1;
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

		if (cp->pJ_LP_PairList == nullptr)
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

PJ_LP proj_helmert_apply(PJ *P, PJ_LP lp, PJ_DIRECTION direction)
{
	struct COMMONPOINTS *cp;
	int inverse;
	PJ_LP out;
	 
	out.lam = HUGE_VAL; out.phi = HUGE_VAL;
	
	cp = find_cp(P->ctx, lp, P->cplist_count, P->cplist);

	if (cp == nullptr)
	{
		if (P->cplist_count == 1 && strcmp(P->cplist[0]->cp_name, "null") == 0)		 
			out = lp;		 	
		else		 
			pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GRID);			
		 	 
		return out;
	}
	    
	int areaId = areaIdPoint(&lp);
	auto closestPoints = findClosestPoints(cp, lp, areaId, direction);

	calculateLscHelmert(P, &lp, &closestPoints, direction);

	out = lp;
	 
	return out;
}

PJ_LP helmert_apply(PJ *P, PJ_LP lp, PJ_DIRECTION direction)
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

PJ_LP collocation_apply(PJ *P, PJ_LP lp, PJ_DIRECTION direction)
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

	if (true)
	{
		struct COMMONPOINTS *cp;
		cp = find_cp(P->ctx, point.lp, P->cplist_count, P->cplist);

		if (cp == nullptr)
		{
			pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GRID);
			return point.xyz;
		}

		int areaId = areaIdPoint(&point.lp);
		auto closestPoints = findClosestPoints(cp, point.lp, areaId, PJ_FWD);
		if (closestPoints.size() == 0)
		{
			pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GRID);
			return point.xyz;
		}		
		if (calculateHelmertParameter(P, &point.lp, &closestPoints, PJ_FWD) == 0);
		{
			pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GRID);
			return point.xyz;
		}
		point.lp = helmert_apply(P, point.lp, PJ_FWD);
		point.lp = collocation_apply(P, point.lp, PJ_FWD);
	}
	else
		point.lp = proj_helmert_apply(P, point.lp, PJ_FWD);

	return point.xyz;
}

static PJ_LPZ reverse_3d(PJ_XYZ xyz, PJ *P)
{
	struct pj_opaque_lschelmert *Q = (struct pj_opaque_lschelmert *) P->opaque;

	PJ_COORD point = { {0,0,0,0} };
	point.xyz = xyz;

	point.lp = proj_helmert_apply(P, point.lp, PJ_INV);

	return point.lpz;
}

PJ *TRANSFORMATION(lschelmert, 0)
{	
	struct pj_opaque_lschelmert *Q = static_cast<struct pj_opaque_lschelmert*>(pj_calloc(1, sizeof(struct pj_opaque_lschelmert)));

	if (Q == nullptr )
	{
		return pj_default_destructor(P, ENOMEM);
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
		proj_log_error(P, "hgridshift: +cp_trans parameter missing.");
		return pj_default_destructor(P, PJD_ERR_NO_ARGS);
	}
	
	proj_cp_init(P, "cp_trans");

	if (proj_errno(P))
	{
		proj_log_error(P, "cp_trans: could not find required cp_tran(s).");
		return pj_default_destructor(P, PJD_ERR_FAILED_TO_LOAD_CPL);
	}
	return P;
}
