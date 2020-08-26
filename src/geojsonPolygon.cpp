/*****************************************************************************
* Project:	PROJ
* Purpose:	GeoJson
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
******************************************************************************/

#ifndef FROM_PROJ_CPP
#define FROM_PROJ_CPP
#endif
#define LRU11_DO_NOT_DEFINE_OUT_OF_CLASS_METHODS
 
#include <fstream>
#include <iostream>

#include "filemanager.hpp"
#include "geojsonPolygon.hpp"
#include "proj/io.hpp"
#include "proj/nn.hpp"
#include "proj/internal/lru_cache.hpp"
#include "proj/internal/include_nlohmann_json.hpp"
#include "proj_json_streaming_writer.hpp"
#include "proj_internal.h"
#include "proj/internal/internal.hpp"

using json = nlohmann::json;
 
NS_PROJ_START

using namespace NS_PROJ::internal;

namespace geoJson
{
	// ---------------------------------------------------------------------------

	static util::BaseObjectNNPtr create(const json & j)
	{
		if (!j.is_object())		 
			throw io::ParsingException("JSON object expected");
		 
		auto type = j["type"];
		if (type == "FeatureCollection")		
			return GeoJsonParser().builtGeoJson(j);
		
		throw io::ParsingException("Unsupported value of \"type\"");
	}

	// ---------------------------------------------------------------------------
	
    util::BaseObjectNNPtr createFromInput(const std::string &text, PJ_CONTEXT *ctx)
	{
		if (!text.empty() && text[0] == '{')
		{
			json j;
			try
			{
				j = json::parse(text);
			}
			catch (const std::exception &e)
			{
				pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_GEOJSON);
				throw io::ParsingException(e.what());
			}
			return create(j);
		}
		throw io::ParsingException("unrecognized format / unknown name");
	}

	// ---------------------------------------------------------------------------	
	
	GeoJson::GeoJson() = default;	 

	// ---------------------------------------------------------------------------
	
	GeoJson::GeoJson(const string &name, const std::vector<FeatureNNPtr> &features, const GeoJsonCrsPtr &crs)		
	{
		m_name = name;
		m_features = features;
		m_GeoJsonCrs = crs;
	};

	// ---------------------------------------------------------------------------

	GeoJson::GeoJson(const string &name, const std::map<std::string, FeatureNNPtr> &features, const GeoJsonCrsPtr &crs)
	{
		m_name = name;
		m_featuresMap = features;
		m_GeoJsonCrs = crs;
	};

	// ---------------------------------------------------------------------------

	GeoJson::~GeoJson() = default;

	// ---------------------------------------------------------------------------

    void GeoJson::_exportToJSON(io::JSONFormatter *formatter) const
	{
		auto objectContext(formatter->MakeObjectContext("FeatureCollection", false));		
		auto writer = formatter->writer();
		
		writer->AddObjKey("name");
		writer->Add(name());

		writer->AddObjKey("crs");
		geoJsonCrs()->exportToJSON(formatter);
		
		writer->AddObjKey("features");

		writer->StartArray();

		for (auto feature : featuresMap())
			feature.second->exportToJSON(formatter);

		writer->EndArray();
		
		// TODO: Delete this
		// Output to console, for testing.	
		std::cout << formatter->toString() << "\n\n";
	}
	
	// ---------------------------------------------------------------------------

	GeoJsonNNPtr GeoJson::create(const string &name, const std::vector<FeatureNNPtr> &features, const GeoJsonCrsPtr &crs)
	{
		auto geojson(GeoJson::nn_make_shared<GeoJson>(name, features, crs));

		return geojson;
	}

	GeoJsonNNPtr GeoJson::create(const string &name, const std::map<std::string, FeatureNNPtr> &features, const GeoJsonCrsPtr &crs)
	{
		auto geojson(GeoJson::nn_make_shared<GeoJson>(name, features, crs));

		return geojson;
	}

	GeoJsonPtr GeoJson::openGeoJson(PJ_CONTEXT *ctx, const std::string &filename)
	{
		auto fp = FileManager::open_resource_file(ctx, filename.c_str());
		
		auto geoJson(GeoJson::nn_make_shared<GeoJson>());
		
		if (!fp)
			return shared_ptr<GeoJson>(nullptr);

		fp->seek(0, SEEK_END);
		unsigned long long fsize = fp->tell();
		fp->seek(0, SEEK_SET);

		char *string = (char *)malloc(fsize + 1);
		fp->read(string, fsize);

		string[fsize] = 0;

		const auto actualName(fp->name());

		if (ends_with(tolower(actualName), "geojson"))
		{
			auto obj = createFromInput(string, ctx);
			auto jsonObj = dropbox::oxygen::nn_dynamic_pointer_cast<GeoJson>(obj);
			
		 	return jsonObj;
		}
		return geoJson;
	}

	// ---------------------------------------------------------------------------

	MultiPolygon::MultiPolygon() = default;

	// ---------------------------------------------------------------------------

	MultiPolygon::MultiPolygon(std::vector<PolygonPoint> &coordinates)
	{
		m_coordinates = coordinates;
	}

	// ---------------------------------------------------------------------------

	MultiPolygon::~MultiPolygon() = default;
	
	// ---------------------------------------------------------------------------
	
	bool MultiPolygon::IsPointInArea(PJ_LP *lp)
	{
		PolygonPoint point = { lp->lam, lp->phi };
		PolygonPoint *vectorPointer = m_coordinates.data();
		int n = (int)size(m_coordinates);

		return isInside(vectorPointer, n, point);
	}

	// ---------------------------------------------------------------------------

	void MultiPolygon::_exportToJSON(io::JSONFormatter *formatter) const
	{
		auto objectContext(formatter->MakeObjectContext("MultiPolygon", false));
		auto writer = formatter->writer();

		writer->AddObjKey("coordinates");

		writer->StartArray();
		writer->StartArray();
		writer->StartArray();

	 	for (auto coordinate : coordinates())
		{
			writer->StartArray();
			writer->Add(newPrecision(coordinate.x, 4), 4);
			writer->Add(newPrecision(coordinate.y, 4), 4);
			writer->EndArray();
		}

		writer->EndArray();
		writer->EndArray();
		writer->EndArray();
	}

	// ---------------------------------------------------------------------------

	MultiPolygonNNPtr MultiPolygon::create(std::vector<PolygonPoint> coordinates)
	{
		auto multiPolygon(MultiPolygon::nn_make_shared<MultiPolygon>(coordinates));

		return multiPolygon;
	}

	// ---------------------------------------------------------------------------

	Feature::Feature() = default;

	// ---------------------------------------------------------------------------

	Feature::Feature(const int &areaid, const MultiPolygonPtr &multipolygonPtr)
	{
		m_s_areaid = to_string(areaid);
		m_areaid = areaid;
		MultiPolygon(multipolygonPtr);
	}

	// ---------------------------------------------------------------------------

	Feature::Feature(const int &areaid, const string &name, const GeoJsonPointPtr &pointPtr)
	{
		m_areaid = areaid;
		m_name = name;
		GeoJsonPoint(pointPtr);
	}

	// ---------------------------------------------------------------------------

	Feature::~Feature() = default;

	// ---------------------------------------------------------------------------
		
	FeatureNNPtr Feature::create(const int &areaid, const MultiPolygonPtr &multipolygonPtr)
	{
		auto feature(Feature::nn_make_shared<Feature>(areaid, multipolygonPtr));	 

		return feature;
	}

	FeatureNNPtr Feature::create(const int &areaid, const string &name, const GeoJsonPointPtr &pointPtr)
	{
		auto feature(Feature::nn_make_shared<Feature>(areaid, name, pointPtr));

		return feature;
	}

	// ---------------------------------------------------------------------------

	void Feature::_exportToJSON(io::JSONFormatter *formatter) const
	{
		auto objectContext(formatter->MakeObjectContext("Feature", false));
		
		auto writer = formatter->writer();

		writer->AddObjKey("properties");
		
		writer->StartObj();
		if (Name() != "")
		{
			writer->AddObjKey("PointName");
			writer->Add(Name());
		}
		if (AreaId() >= 0)
		{
			writer->AddObjKey("areaid");
			writer->Add(AreaId());
		}
		writer->EndObj();

		writer->AddObjKey("geometry");

		if (MultiPolygon() != nullptr)
			MultiPolygon()->exportToJSON(formatter);
		else if (Point() != nullptr)
			Point()->exportToJSON(formatter);
	}

	// ---------------------------------------------------------------------------
	
	GeoJsonPoint::GeoJsonPoint(const double &x, const double &y)
	{
		m_x = x;
		m_y = y;
	}

	// ---------------------------------------------------------------------------

	GeoJsonPoint::~GeoJsonPoint() = default;

	// ---------------------------------------------------------------------------

	void GeoJsonPoint::_exportToJSON(io::JSONFormatter *formatter) const
	{
		auto objectContext(formatter->MakeObjectContext("Point", false));
		auto writer = formatter->writer();

		writer->AddObjKey("coordinates");
	
		writer->StartArray();
		writer->Add(X(), 14);
		writer->Add(Y(), 14);
		writer->EndArray();
	}
	
	// ---------------------------------------------------------------------------

	GeoJsonPointNNPtr GeoJsonPoint::create(const double &x, const double &y)
	{
		GeoJsonPointNNPtr geojsonpoint(GeoJsonPoint::nn_make_shared<GeoJsonPoint>(x, y));

		return geojsonpoint;
	}

	// ---------------------------------------------------------------------------
	
	GeoJsonCrs::GeoJsonCrs(const string &name)
	{
		m_name = name;
	}
	
	// ---------------------------------------------------------------------------
	 
	GeoJsonCrs::~GeoJsonCrs() = default;

	// ---------------------------------------------------------------------------

	void GeoJsonCrs::_exportToJSON(io::JSONFormatter *formatter) const
	{
		auto objectContext(formatter->MakeObjectContext("name", false));		
		auto writer = formatter->writer();

		writer->AddObjKey("properties");

		writer->StartObj();
		writer->AddObjKey("name");
		writer->Add(Name());
		writer->EndObj();
	}

	// ---------------------------------------------------------------------------

	GeoJsonCrsNNPtr GeoJsonCrs::create(const string &name)
	{
		auto geojsoncrs(GeoJsonCrs::nn_make_shared<GeoJsonCrs>(name));

		return geojsoncrs;
	}

	// ---------------------------------------------------------------------------	

	GeoJsonGeometry::GeoJsonGeometry() = default;

	// ---------------------------------------------------------------------------

	GeoJsonGeometry::~GeoJsonGeometry() = default;
	   
	// ---------------------------------------------------------------------------
	
	GeoJsonNNPtr GeoJsonParser::builtGeoJson(const json &j)
	{
		GeoJsonPtr geoJson;
		 
		auto type = j["type"].get<std::string>();
		
		if (type == "FeatureCollection")
		{
			auto featureCollection = j["type"];

			if (!j.contains("name"))
				throw io::ParsingException("\"name\" is missing");
			
			std::string name = j["name"].get<std::string>();

			if (!j.contains("crs"))
				throw io::ParsingException("\"crs\" is missing");			 
			 
			std::map<std::string, FeatureNNPtr> featurePtrMap;

			auto crs(GeoJsonCrs::nn_make_shared<GeoJsonCrs>(""));

			if (j.contains("crs"))
			{
				auto crsJson = j["crs"];

				crs = builtCrs(crsJson);
			}
	 		if (j.contains("features"))
			{
				auto features = j["features"];				 

				for (const auto &feature : features)
				{
					if (feature["type"].get<std::string>() == "Feature")
					{
						auto feat = builtFeature(feature);

						if (feat->Name() != "")
							featurePtrMap.emplace(feat->Name(), feat);
						else if (feat->AreaIdString() != "")
							featurePtrMap.emplace(feat->AreaIdString(), feat);
					}
				}
			}
			return GeoJson::create(name, featurePtrMap, crs);
		}
		throw io::ParsingException("builtGeoJson failed.");
	} 
 
	// ---------------------------------------------------------------------------
	 
	MultiPolygonNNPtr GeoJsonParser::builtMultipolygon(const json &j)
	{
		std::vector<PolygonPoint> coordinateList{};

		if (!j.contains("coordinates"))
			throw io::ParsingException("builtMultipolygon failed.");

		auto coordinates = j["coordinates"];

		for (const auto &coordinate : coordinates)			 
			recursive_iterate(coordinate, coordinateList, [](json::const_iterator) {});
		
		return MultiPolygon::create(coordinateList);
	}
	 
	// ---------------------------------------------------------------------------
	 
	GeoJsonPointNNPtr GeoJsonParser::builtPoint(const json &j)
	{
		if (!j.contains("coordinates"))
			throw io::ParsingException("builtPoint failed.");

		auto coordinates = j["coordinates"];

		double x = 0.0;
		double y = 0.0;

		if (coordinates.is_array())
		{
			for (auto it = coordinates.begin(); it != coordinates.end(); ++it)
			{
				if (it->is_number_float())
				{
					x = it.value();

					if (it != coordinates.end())
					{
						++it;

						if (it->is_number_float())
							y = it.value();
					}
					break;
				}
			}
		}
		return GeoJsonPoint::create(x, y);
	}
	 
	// ---------------------------------------------------------------------------
	 
	FeatureNNPtr GeoJsonParser::builtFeature(const json &j)
	{
		MultiPolygonPtr multipolygonPtr;
		GeoJsonPointPtr geoJsonPointPtr;

		if (j.contains("properties"))
		{
			auto properties = j["properties"];

			if (!properties.contains("areaid"))
				throw io::ParsingException("\"areaid\" is missing");

			int areaid = (int)properties["areaid"].get<double>();
		  
			if (!j.contains("geometry"))
				throw io::ParsingException("\"geometry\" is missing");

			auto geometry = j["geometry"];
			 
			if (geometry["type"].get<std::string>() == "MultiPolygon")
			{
				multipolygonPtr = util::nn_dynamic_pointer_cast<MultiPolygon>(builtMultipolygon(geometry));
				return Feature::create(areaid, multipolygonPtr);
			}
			else if (geometry["type"].get<std::string>() == "Point")
			{
				if (!properties.contains("PointName"))
					throw io::ParsingException("Point name is missing.");
				
				auto pointName = properties["PointName"].get<std::string>();

				geoJsonPointPtr = util::nn_dynamic_pointer_cast<GeoJsonPoint>(builtPoint(geometry));
				return Feature::create(areaid, pointName, geoJsonPointPtr);
			}
			throw io::ParsingException("built method missing.");
		}
		throw io::ParsingException("builtFeature failed.");
	}

	// ---------------------------------------------------------------------------

	GeoJsonCrsNNPtr GeoJsonParser::builtCrs(const json &j)
	{
		if (!j.contains("properties"))
			throw io::ParsingException("builtCrs failed.");

		auto crs(GeoJsonCrs::nn_make_shared<GeoJsonCrs>(""));

		auto prop = j["properties"];

		if (prop.contains("name"))
		{
			auto name = prop["name"].get<std::string>();		 

			return GeoJsonCrs::create(name);
		}
		return crs;
	}

	util::BaseObjectNNPtr GeoJsonParser::create(const json & j)
	{
		if (!j.is_object())
			throw io::ParsingException("JSON object expected");

		auto type = j["type"];
		if (type == "FeatureCollection")
			return GeoJsonParser().builtGeoJson(j);

		throw io::ParsingException("Unsupported value of \"type\"");
	}	 
	 
	// ---------------------------------------------------------------------------

	ListOfGeoJson pj_geojson_init(PJ *P, const char *geojsonkey)
	{
		std::string key("s");
		key += geojsonkey;

		const char *geojsonnames = pj_param(P->ctx, P->params, key.c_str()).s;

		if (geojsonnames == nullptr)
			return {};

		auto listOfGeoJsonNames = split(std::string(geojsonnames), ',');

		ListOfGeoJson geoJsons;

		for (const auto &geoJsonStr : listOfGeoJsonNames)
		{
			const char *geoJsonname = geoJsonStr.c_str();

			if (geoJsonname[0] == '@')
				geoJsonname++;
			 
			auto polyset = geoJson::GeoJson::openGeoJson(P->ctx, geoJsonname);

			if (!polyset)
			{
				pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GEOJSON);

				return {};
			}
			else
				geoJsons.emplace_back(std::move(polyset));
		}
		return geoJsons;
	}

	// ---------------------------------------------------------------------------
	 
	int areaIdPoint(PJ *P, const ListOfGeoJson &geoJsonList, PJ_LP *lp)
	{
		for (const auto& geoJson : geoJsonList)
		{
			for (const auto& feature : geoJson->featuresMap())
			{
				if (!feature.second.get())
				{
					if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
						proj_log_trace(P, "GeoJson feature is empty");

					return 0;
				}
				
				auto feat = feature.second;
				auto areaId = feat->AreaId();
			 
				if (!feat->MultiPolygon())
				{
					if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
						proj_log_trace(P, "GeoJson feature is empty");

					return 0;
				}
				if (feat->MultiPolygon()->IsPointInArea(lp))
				{
					if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
						proj_log_trace(P, "Input point was found in area with id %d", areaId);

					return areaId;
				}
			}
		}
		if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE)
			proj_log_trace(P, "Input point was not found in any areas");

		return 0;
	}
}
NS_PROJ_END