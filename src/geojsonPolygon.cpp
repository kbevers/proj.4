/*****************************************************************************
* Project:	PROJ
* Purpose:	GeoJson Multipolygon
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

#ifndef FROM_PROJ_CPP
#define FROM_PROJ_CPP
#endif
#define LRU11_DO_NOT_DEFINE_OUT_OF_CLASS_METHODS
 
#include <fstream> // std::ifstream
#include <iostream>

#include "filemanager.hpp"
#include "geojsonPolygon.hpp"
#include "proj/internal/lru_cache.hpp"
#include "proj/internal/nlohmann/json.hpp"
#include "proj_internal.h"
#include "proj/internal/internal.hpp" // for split
 
using json = nlohmann::json;
 
NS_PROJ_START

using namespace NS_PROJ::internal;

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

					// Converts to radians
					PolygonPoint p{ proj_torad(x), proj_torad(y) };
					vlist.push_back(p);
				}
			}
		}
		else
			f(it);
	}
}

std::unique_ptr<GeoJsonMultiPolygonSet>
GeoJsonMultiPolygonSet::open(PJ_CONTEXT *ctx, const std::string &filename)
{
	if (filename == "null")
	{
		auto polygonSet = std::unique_ptr<GeoJsonMultiPolygonSet>(new GeoJsonMultiPolygonSet());
		polygonSet->m_name = filename;
		polygonSet->m_format = "null";
		//set->m_polygons.push_back(std::unique_ptr<NullVerticalShiftGrid>(new NullVerticalShiftGrid()));
		return polygonSet;
	}
	auto fp = FileManager::open_resource_file(ctx, filename.c_str());
	if (!fp)
		return nullptr;
	 
	const auto actualName(fp->name());
	
	if (ends_with(tolower(actualName), "geojson"))
	{	 
		auto polygonSet = GeoJsonMultiPolygonSet::parse(ctx, std::move(fp), actualName);
		//auto polygon = GeoJsonMultiPolygon::open(ctx, std::move(fp), actualName);

	 	if (!polygonSet)
	 		return nullptr;

		//auto set = std::unique_ptr<GeoJsonMultiPolygonSet>(new GeoJsonMultiPolygonSet());	
		polygonSet->m_format = "geojson";
		polygonSet->m_name = filename;

		return polygonSet;
	}
};

bool GeoJsonMultiPolygonSet::reopen(PJ_CONTEXT *ctx) 
{
	pj_log(ctx, PJ_LOG_DEBUG_MAJOR, "Polygon %s has changed. Re-loading it", m_name.c_str());

	auto newGS = open(ctx, m_name);
	m_format.clear();
	
	if (newGS)
	{
		// TODO: Complete this...
		//	m_grids = std::move(newGS->m_grids);
	}
	return !m_format.empty();
}

std::unique_ptr<GeoJsonMultiPolygonSet>
GeoJsonMultiPolygonSet::parse(PJ_CONTEXT *ctx, std::unique_ptr<File> fp, const std::string &filename)
{
	auto set = std::unique_ptr<GeoJsonMultiPolygonSet>(new GeoJsonMultiPolygonSet());

	fp->seek(0, SEEK_END);
	long fsize = fp->tell();
	fp->seek(0, SEEK_SET);

	char *string = (char *)malloc(fsize + 1);
	fp->read(string, fsize);
	
	string[fsize] = 0;

	// parse and serialize JSON
	json j_complete = json::parse(string);
	
	// Output to console. For testing.
	// std::cout << std::setw(4) << j_complete << "\n\n";

	auto feat = j_complete.at("features");

	for (auto it = feat.begin(); it != feat.end(); ++it)
	{	
		bool isMultiPolygon = false;
		vector<PolygonPoint> pointVector;

		auto areas = (*it)["properties"];
		auto area = areas.find("areaid"); 
		
		if (area.value().is_number_integer())
		{
			__int32 id = area.value();
			auto polygon = new GeoJsonMultiPolygon(id);
			
			auto geo = (*it)["geometry"];
			
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
					polygon->m_pointList = pointVector;
				}
				if (isMultiPolygon)
				{
					// set->m_polygons->push_back(polygon);
					set->m_polygons.push_back(std::unique_ptr<GeoJsonMultiPolygon>(polygon));
				}
			}		
		}
	}

	// Testing Json dump
	auto json_string = j_complete.dump();
	
	return set;
}

/*
void GeoJsonMultiPolygonSet::reassign_context(PJ_CONTEXT *ctx)
{
	for (const auto &poly : m_polygons)
	{
		poly->reassign_context(ctx);
	}
} 
*/

GeoJsonMultiPolygonSet::GeoJsonMultiPolygonSet() = default;

GeoJsonMultiPolygonSet::~GeoJsonMultiPolygonSet() = default;

Polygon::Polygon(const __int32 &areaid) {};

Polygon::~Polygon() = default;

GeoJsonMultiPolygon::GeoJsonMultiPolygon(__int32 &areaid) : Polygon(areaid)
{
	m_areaid = areaid;
};

GeoJsonMultiPolygon::~GeoJsonMultiPolygon() = default;

GeoJsonMultiPolygon *GeoJsonMultiPolygon::open(PJ_CONTEXT *ctx, std::unique_ptr<File> fp, const std::string &name)
{
	const char *cstr = name.c_str();

	FILE *f = fopen(cstr, "rb");
	 
	__int32 testId = 2;
	auto set = new GeoJsonMultiPolygon(testId/*name*/);

	fclose(f);

	return set;
}
 
ListOfMultiPolygons pj_polygon_init(PJ *P, const char *polygonkey)
{
	std::string key("s");
	key += polygonkey;
	const char *polygonnames = pj_param(P->ctx, P->params, key.c_str()).s;

	if (polygonnames == nullptr)
		return {};

	auto listOfPolygonNames = split(std::string(polygonnames), ',');
	ListOfMultiPolygons polygons;

	for (const auto &polygonStr : listOfPolygonNames)
	{
		const char *polygonname = polygonStr.c_str();
		bool canFail = false;
		if (polygonname[0] == '@')
		{
			canFail = true;
			polygonname++;
		}	
	    auto polySet = GeoJsonMultiPolygonSet::open(P->ctx, polygonname);
		
		if (!polySet)
		{
	 		if (!canFail) 
			{
			   // if (proj_context_errno(P->ctx) != PJD_ERR_NETWORK_ERROR)
				{
					pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GEOJSON);
				}
				return {};
			}
			pj_ctx_set_errno(P->ctx, 0);
		}
		else 
		{
			polygons.emplace_back(std::move(polySet)); 
		}
	}
	return polygons;
}

bool pointIsInAreaTest(PJ_LP *lp, GeoJsonMultiPolygon &polygon)
{
	PolygonPoint point = { lp->lam, lp->phi}; // TODO: Feil eining

	PolygonPoint *vectorPointer = polygon.m_pointList.data();
	
	int n = (int)size(polygon.m_pointList);
	
	return isInside(vectorPointer, n, point);	
}

bool pointIsInArea(PJ_LP pointPJ_LP, char* fileName)
{
	std::ifstream file(fileName, std::ios::in);

	PolygonPoint points[] = { {} };
	PolygonPoint point = { pointPJ_LP.phi, pointPJ_LP.lam }; // TODO: Feil eining
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

__int32 areaIdPoint(const ListOfMultiPolygons &polygonList, PJ_LP *lp)
{
	for (const auto& polygonSet : polygonList)
	{
		// TODO: This is dirty. Refactorize.
		for (auto polygon = polygonSet->polygons().begin(); polygon != polygonSet->polygons().end(); polygon++)
		{			
			if (static_cast<GeoJsonMultiPolygon *>(polygon->get()))
			{
				auto poly = polygon->get();
				int id = poly->Id();

				if (pointIsInAreaTest(lp, *poly))
				   return id;
			}
			//std::cout << *polygon << std::endl;
		}
	} 
	return 0;
}
NS_PROJ_END