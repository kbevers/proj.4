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
******************************************************************************/

#include <vector>
#include <iostream>

#include "proj.h"
#include "proj/common.hpp"
#include "proj/io.hpp"
#include "proj/nn.hpp"
#include "proj/util.hpp"
#include "proj/internal/include_nlohmann_json.hpp"
#include "point_in_polygon.h"
#include "filemanager.hpp"

using json = nlohmann::json;

NS_PROJ_START

// ---------------------------------------------------------------------------

namespace geoJson
{ 
	PROJ_DLL util::BaseObjectNNPtr createFromInput(const std::string &text, PJ_CONTEXT *ctx);

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
				double x = it.value();
				if (it != j.end())
				{
					++it;
					if (it->is_number_float())
					{
						double y = it.value();

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
	 
	// ---------------------------------------------------------------------------
	 
	inline double newPrecision(double n, int i)
	{
		double num = (float)(round(pow(10, i) * n));
		int denom = (int)(pow(10, i));

		return num / denom;
	}
 
	class GeoJsonCrs;
	using GeoJsonCrsPtr = std::shared_ptr<GeoJsonCrs>;
	using GeoJsonCrsNNPtr = util::nn<GeoJsonCrsPtr>;

	class PROJ_GCC_DLL GeoJsonCrs final :
	    public util::BaseObject,
		public io::IJSONExportable
	{
	protected:
		string m_name{};
	public:		
		PROJ_FOR_TEST explicit GeoJsonCrs(const string &name);
		INLINED_MAKE_SHARED

	    PROJ_INTERNAL ~GeoJsonCrs() override;
		PROJ_FOR_TEST const std::string &Name() const { return m_name; }
		PROJ_FOR_TEST void Name(std::string x) { m_name = std::move(x); }
		PROJ_INTERNAL void _exportToJSON(io::JSONFormatter *formatter) const override;
		PROJ_INTERNAL static GeoJsonCrsNNPtr create(const string &name);
	};

	// ---------------------------------------------------------------------------

	class GeoJsonGeometry;
	using GeoJsonGeometryPtr = std::shared_ptr<GeoJsonGeometry>;
	using GeoJsonGeometryNNPtr = util::nn<GeoJsonGeometryPtr>;

	class PROJ_GCC_DLL GeoJsonGeometry :
		public util::BaseObject,
		public io::IJSONExportable
	{
	protected:
		PROJ_INTERNAL GeoJsonGeometry();
	public:
		PROJ_INTERNAL ~GeoJsonGeometry() override;
	}; 

	// ---------------------------------------------------------------------------

	class GeoJsonPoint;
	using GeoJsonPointPtr = std::shared_ptr<GeoJsonPoint>;
	using GeoJsonPointNNPtr = util::nn<GeoJsonPointPtr>;
 
	class PROJ_GCC_DLL GeoJsonPoint final :
		public GeoJsonGeometry	   
	{
	protected:
		double m_x = 0;
		double m_y = 0;
	public:
		PROJ_FOR_TEST explicit GeoJsonPoint(const double &x, const double &y);
		INLINED_MAKE_SHARED

		PROJ_FOR_TEST ~GeoJsonPoint() override;

		PROJ_FOR_TEST double X() const { return m_x; }
		PROJ_FOR_TEST void X(double x) { m_x = std::move(x); }

		PROJ_FOR_TEST double Y() const { return m_y; }
		PROJ_FOR_TEST void Y(double y) { m_y = std::move(y); }

		PROJ_FOR_TEST double X_rad() const { return proj_torad(m_x); }
		PROJ_FOR_TEST double Y_rad() const { return proj_torad(m_y); }

		PROJ_INTERNAL void _exportToJSON(io::JSONFormatter *formatter) const override;
		PROJ_INTERNAL static GeoJsonPointNNPtr create(const double &x, const double &y);
	};
	
	// ---------------------------------------------------------------------------

	class MultiPolygon;
	using MultiPolygonPtr = std::shared_ptr<MultiPolygon>;
	using MultiPolygonNNPtr = util::nn<MultiPolygonPtr>;

	class PROJ_GCC_DLL MultiPolygon final :
	    public GeoJsonGeometry
	{
	protected:
		std::vector<PolygonPoint> m_coordinates{};
	public:
		PROJ_FOR_TEST MultiPolygon();

		PROJ_FOR_TEST explicit MultiPolygon(std::vector<PolygonPoint> &coordinates);
		INLINED_MAKE_SHARED

		PROJ_FOR_TEST ~MultiPolygon() override;

		PROJ_FOR_TEST const std::vector<PolygonPoint> &coordinates() const { return m_coordinates; }
		PROJ_FOR_TEST bool IsPointInArea(PJ_LP *lp);
		PROJ_INTERNAL void _exportToJSON(io::JSONFormatter *formatter) const override;
		PROJ_INTERNAL static MultiPolygonNNPtr create(std::vector<PolygonPoint> coordinates);
	};

	// ---------------------------------------------------------------------------
	
	class Feature;
	using FeaturePtr = std::shared_ptr<Feature>;
	using FeatureNNPtr = util::nn<FeaturePtr>;

	class PROJ_GCC_DLL Feature final :
		public util::BaseObject,
		public io::IJSONExportable
	{
	protected:
		int m_areaid = -1;
		string m_name{};
		string m_s_areaid{};
		MultiPolygonPtr m_multipolygonPtr;
		GeoJsonPointPtr m_pointPtr;
		std::vector<std::unique_ptr<PJ_LPZ>> m_coordinates{};
	public:
		PROJ_FOR_TEST Feature();

		PROJ_FOR_TEST explicit Feature(const int &areaid, const MultiPolygonPtr &multipolygonPtr);
		PROJ_FOR_TEST explicit Feature(const int &areaid, const string &name, const GeoJsonPointPtr &pointPtr);
		INLINED_MAKE_SHARED
		
		PROJ_FOR_TEST ~Feature() override;

		PROJ_FOR_TEST const int &AreaId() const { return m_areaid; }
		PROJ_FOR_TEST const std::string &Name() const { return m_name; }
		PROJ_FOR_TEST const std::string &AreaIdString() const { return m_s_areaid; }
		PROJ_FOR_TEST void Name(std::string x) { m_name = std::move(x); }
		
		PROJ_FOR_TEST const MultiPolygonPtr &MultiPolygon() const { return m_multipolygonPtr; };
		PROJ_FOR_TEST void MultiPolygon(MultiPolygonPtr multipolygonPtr) { m_multipolygonPtr = std::move(multipolygonPtr); }
		PROJ_FOR_TEST const GeoJsonPointPtr &Point() const { return m_pointPtr; };
		PROJ_FOR_TEST void GeoJsonPoint(GeoJsonPointPtr pointPtr) { m_pointPtr = std::move(pointPtr); }

		PROJ_INTERNAL void _exportToJSON(io::JSONFormatter *formatter) const override;
		PROJ_INTERNAL static FeatureNNPtr create(const int &areaid, const MultiPolygonPtr &multipolygonPtr);
		PROJ_INTERNAL static FeatureNNPtr create(const int &areaid, const string &name, const GeoJsonPointPtr &pointPtr);
		
		PROJ_FOR_TEST bool operator==(const Feature & obj) const
		{
			 if (Name() == obj.Name())
				return true;

			return false;
		};		 	
	};

	// ---------------------------------------------------------------------------

	class GeoJson;
	using GeoJsonPtr = std::shared_ptr<GeoJson>;
	using GeoJsonNNPtr = util::nn<GeoJsonPtr>;
	
	class PROJ_GCC_DLL GeoJson final :
	public util::BaseObject,
		public io::IJSONExportable
	{
	protected:
		string m_name{};
		GeoJsonCrsPtr m_GeoJsonCrs{};
		std::vector<FeatureNNPtr> m_features{};
		std::map<std::string, FeatureNNPtr> m_featuresMap{};
	public:
		PROJ_INTERNAL GeoJson();
		
		PROJ_FOR_TEST explicit GeoJson(const string &name, const std::vector<FeatureNNPtr> &features, const GeoJsonCrsPtr &crs);
		PROJ_FOR_TEST explicit GeoJson(const string &name, const std::map<std::string, FeatureNNPtr> &features, const GeoJsonCrsPtr &crs);
		INLINED_MAKE_SHARED

		PROJ_INTERNAL ~GeoJson() override;

		PROJ_FOR_TEST const std::string &name() const { return m_name; }
		PROJ_FOR_TEST const std::vector<FeatureNNPtr> &features() const { return m_features; }
		PROJ_FOR_TEST const std::map<std::string, FeatureNNPtr> &featuresMap() const { return m_featuresMap; }
		PROJ_FOR_TEST const GeoJsonCrsPtr &geoJsonCrs() const { return m_GeoJsonCrs; }
		PROJ_INTERNAL void _exportToJSON(io::JSONFormatter *formatter) const override;
		PROJ_INTERNAL static GeoJsonNNPtr create(const string &name, const std::vector<FeatureNNPtr> &features, const GeoJsonCrsPtr &crs);
		PROJ_INTERNAL static GeoJsonNNPtr create(const string &name, const std::map<std::string, FeatureNNPtr> &features, const GeoJsonCrsPtr &crs);
		PROJ_FOR_TEST static GeoJsonPtr openGeoJson(PJ_CONTEXT *ctx, const std::string &filename);
		
		PROJ_FOR_TEST bool FeatureExits(std::string name)
		{
			std::map<std::string, FeatureNNPtr>::iterator it;
			it = m_featuresMap.find(name);
			
			return it != m_featuresMap.end(); 
		};
		PROJ_FOR_TEST const FeatureNNPtr &GetFeature(std::string name)
		{ 
			std::map<std::string, FeatureNNPtr>::iterator it;
			it = m_featuresMap.find(name);
			
			return it->second;		
		};
	};
		 	
	// ---------------------------------------------------------------------------

	class GeoJsonParser
	{
	protected:
	public:
		GeoJsonParser() = default;

		GeoJsonNNPtr builtGeoJson(const json &j);
		MultiPolygonNNPtr builtMultipolygon(const json &j);
		GeoJsonPointNNPtr builtPoint(const json &j);
		FeatureNNPtr builtFeature(const json &j);
		GeoJsonCrsNNPtr builtCrs(const json &j);

		util::BaseObjectNNPtr create(const json &j);
	};

	// ---------------------------------------------------------------------------
	 
	typedef std::vector<std::shared_ptr<GeoJson>> ListOfGeoJson;

	// ---------------------------------------------------------------------------
	 
	ListOfGeoJson pj_geojson_init(PJ *P, const char *polygonkey);

	// ---------------------------------------------------------------------------

	int areaIdPoint(PJ *P, const ListOfGeoJson &geoJsonList, PJ_LP *lp);
}
NS_PROJ_END