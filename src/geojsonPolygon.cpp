/*****************************************************************************
* Project:	PROJ
* Purpose:	GeoJson Multipolygon
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

#ifndef FROM_PROJ_CPP
#define FROM_PROJ_CPP
#endif
#define LRU11_DO_NOT_DEFINE_OUT_OF_CLASS_METHODS
 
#include <fstream> // std::ifstream
#include <iostream>

#include "geojsonPolygon.hpp"
#include "proj/internal/lru_cache.hpp"
#include "proj/internal/nlohmann/json.hpp"
#include "proj_internal.h"
#include "proj/internal/internal.hpp" // for split
 
using json = nlohmann::json;
 
NS_PROJ_START

using namespace NS_PROJ::internal;
 
std::unique_ptr<GeoJsonMultiPolygonSet>
GeoJsonMultiPolygonSet::open(PJ_CONTEXT *ctx, const std::string &filename)
{
	auto set = std::unique_ptr<GeoJsonMultiPolygonSet>(new GeoJsonMultiPolygonSet());
	set->m_name = filename;
	set->m_format = "null";
 
	return set;
};

bool GeoJsonMultiPolygonSet::reopen(PJ_CONTEXT *ctx) 
{
	pj_log(ctx, PJ_LOG_DEBUG_MAJOR, "Grid %s has changed. Re-loading it",
		m_name.c_str());
	auto newGS = open(ctx, m_name);
	m_format.clear();
	if (newGS) {
	//	m_grids = std::move(newGS->m_grids);
	}
	return !m_format.empty();
}

void GeoJsonMultiPolygonSet::reassign_context(PJ_CONTEXT *ctx) {
	for (const auto &poly : m_polygons)
	{
		poly->reassign_context(ctx);
	}
} 

GeoJsonMultiPolygonSet::GeoJsonMultiPolygonSet() = default;

GeoJsonMultiPolygonSet::~GeoJsonMultiPolygonSet() = default;

ListOfMultiPolygon pj_polygon_init(PJ *P, const char *polygonkey)
{
	std::string key("s");
	key += polygonkey;
	const char *polygonnames = pj_param(P->ctx, P->params, key.c_str()).s;

	if (polygonnames == nullptr)
		return {};

	auto listOfPolygonNames = split(std::string(polygonnames), ',');
	ListOfMultiPolygon polygons;

	for (const auto &polygonStr : listOfPolygonNames)
	{
		const char *polygonname = polygonStr.c_str();
		bool canFail = false;
		if (polygonname[0] == '@')
		{
			canFail = true;
			polygonname++;
		}
	
	 /*   auto polySet = GeoJsonMultiPolygonSet::open(P->ctx, polygonname);
			if (!polySet)
		{
			if (!canFail) 
			{
			    if (proj_context_errno(P->ctx) != PJD_ERR_NETWORK_ERROR)
				{
					pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GRID);
				}
				return {};
			}
			pj_ctx_set_errno(P->ctx, 0);
		}
		else 
		{
			polygons.emplace_back(std::move(polySet));
		} */
	} 
	return polygons;
}
NS_PROJ_END

bool pointIsInArea(PJ_LP pointPJ_LP, char* fileName)
{
	std::ifstream file(fileName, std::ios::in);

	PolygonPoint points[] = { {} };
	PolygonPoint point = { pointPJ_LP.phi, pointPJ_LP.lam };
	vector<PolygonPoint> pointVector;

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

			PolygonPoint areaPoint;
			areaPoint.x = x;
			areaPoint.y = y;

			pointVector.push_back(areaPoint);
		}
		file.close();
	};

	PolygonPoint *vectorPointer = pointVector.data();
	int n = size(pointVector);

	return isInside(vectorPointer, n, point);
}

template<class UnaryFunction>
void recursive_iterate(const json& j, vector<PolygonPoint> &vlist, UnaryFunction f)
{
	for (auto it = j.begin(); it != j.end(); ++it)
	{
		auto v = it.value();

		if (it->is_array() || it->is_object())
			recursive_iterate(*it, vlist, f);
		else if (it->is_null())
			f(it);
		else if (it->is_number_float())
		{
			float x = it.value();
			if (it != j.end())
			{
				++it;
				if (it->is_number_float())
				{
					float y = it.value();
					PolygonPoint p{ x,  y };
					vlist.push_back(p);
				}
			}
		}
		else
			f(it);
	}
}
 
int areaIdPoint(PJ_LP *lp)
{
	// TODO: Flytte områdefilene
	char* fileName2 = "C:/Prosjekter/SkTrans/EurefNgo/Punksky_tilfeldig/Area2.csv";
	char* fileName3 = "C:/Prosjekter/SkTrans/EurefNgo/Punksky_tilfeldig/Area3.csv";
	char* fileName4 = "C:/Prosjekter/SkTrans/EurefNgo/Punksky_tilfeldig/Area4.csv";

	// TODO: Områdefil som geojson. Json ligg under include/proj/internal/nlohmann
	// testReadGeojson();

	if (pointIsInArea(*lp, fileName2))
		return 2;
	else if (pointIsInArea(*lp, fileName3))
		return 3;
	else if (pointIsInArea(*lp, fileName4))
		return 4;

	return 1;
}

void testReadGeojson(/*char* fileName*/)
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
	//std::cout << std::setw(4) << j_complete << "\n\n";

	//recursive_iterate(j_complete, [](json::const_iterator it)
	//{});

	// Test
	auto feat = j_complete.at("features");

	for (auto it1 = feat.begin(); it1 != feat.end(); ++it1)
	{
		bool isMultiPolygon = false;
		vector<PolygonPoint> pointVector;

		auto geo = (*it1)["geometry"];

		for (auto& el : geo.items())
		{
			if (el.key() == "type")
			{
				if (el.value() == "MultiPolygon")
					isMultiPolygon = true;
			}
			if (el.key() == "coordinates")
			{
				recursive_iterate(el, pointVector, [](json::const_iterator it) {});
			}
			if (isMultiPolygon)
			{
			}
		}
	}
	auto json_string = j_complete.dump();
}
