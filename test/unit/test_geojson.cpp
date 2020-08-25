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

#include "gtest_include.h"

#define FROM_PROJ_CPP

#include "proj/common.hpp"
#include "proj/datum.hpp"
#include "proj/io.hpp"
#include "proj/metadata.hpp"
#include "proj/util.hpp"
#include "..\src\geojsonPolygon.hpp"

#include "proj/internal/internal.hpp"

using namespace osgeo::proj::common;
using namespace osgeo::proj::datum;
using namespace osgeo::proj::io;
using namespace osgeo::proj::metadata;
using namespace osgeo::proj::util;
using namespace osgeo::proj::geoJson;

TEST(json_import, geojsonmultipolygon)
{
	auto json =
		"{\n"
		"  \"$schema\": \"foo\",\n"
		"  \"type\": \"FeatureCollection\",\n"
		"  \"name\": \"NGO_areas\",\n"
		"  \"crs\": {\n"
		"    \"type\": \"name\",\n"
		"    \"properties\": {\n"
		"      \"name\": \"urn:ogc:def:crs:OGC:1.3:CRS84\"\n"
		"    }\n"
		"  },\n"
		"  \"features\": [\n"
		"    {\n"
		"      \"type\": \"Feature\",\n"
		"      \"properties\": {\n"
		"        \"areaid\": 1\n"
		"      },\n"
		"      \"geometry\": {\n"
		"        \"type\": \"MultiPolygon\",\n"
		"        \"coordinates\": [\n"
		"          [\n"
		"            [\n"
		"              [\n"
		"                9.17,\n"
		"                58.87\n"
		"              ],\n"
		"              [\n"
		"                9.17,\n"
		"                58.88\n"
		"              ],\n"
		"              [\n"
		"                9.18,\n"
		"                58.88\n"
		"              ],\n"
		"              [\n"
		"                9.18,\n"
		"                58.87\n"
		"              ]\n"
		"            ]\n"
		"          ]\n"
		"        ]\n"
		"      }\n"
		"    },\n"
		"    {\n"
		"      \"type\": \"Feature\",\n"
		"      \"properties\": {\n"
		"        \"areaid\": 2\n"
		"      },\n"
		"      \"geometry\": {\n"
		"        \"type\": \"MultiPolygon\",\n"
		"        \"coordinates\": [\n"
		"          [\n"
		"            [\n"
		"              [\n"
		"                9.19,\n"
		"                58.85\n"
		"              ],\n"
		"              [\n"
		"                9.19,\n"
		"                58.86\n"
		"              ],\n"
		"              [\n"
		"                9.21,\n"
		"                58.86\n"
		"              ],\n"
		"              [\n"
		"                9.21,\n"
		"                58.85\n"
		"              ]\n"
		"            ]\n"
		"          ]\n"
		"        ]\n"
		"      }\n"
		"    }\n"
		"  ]\n"
		"}";
 
	auto obj = createFromInput(json, nullptr);
	ASSERT_TRUE(obj != nullptr);

	auto jsonObj = nn_dynamic_pointer_cast<GeoJson>(obj);
	ASSERT_TRUE(jsonObj != nullptr);

	EXPECT_EQ(jsonObj->exportToJSON(&(JSONFormatter::create()->setSchema("foo"))), json);
}

TEST(json_import, geojsonpoints)
{
	auto json =
		"{\n"
		"  \"$schema\": \"foo\",\n"
		"  \"type\": \"FeatureCollection\",\n"
		"  \"name\": \"SourcePoints\",\n"
		"  \"crs\": {\n"
		"    \"type\": \"name\",\n"
		"    \"properties\": {\n"
		"      \"name\": \"urn:ogc:def:crs:EPSG::4258\"\n"
		"    }\n"
		"  },\n"
		"  \"features\": [\n"
		"    {\n"
		"      \"type\": \"Feature\",\n"
		"      \"properties\": {\n"
		"        \"PointName\": \"D41T0006\",\n"
		"        \"areaid\": 1\n"
		"      },\n"
		"      \"geometry\": {\n"
		"        \"type\": \"Point\",\n"
		"        \"coordinates\": [\n"
		"          7.500665,\n"
		"          57.967922\n"
		"        ]\n"
		"      }\n"
		"    },\n"
		"    {\n"
		"      \"type\": \"Feature\",\n"
		"      \"properties\": {\n"
		"        \"PointName\": \"D41T0015\",\n"
		"        \"areaid\": 1\n"
		"      },\n"
		"      \"geometry\": {\n"
		"        \"type\": \"Point\",\n"
		"        \"coordinates\": [\n"
		"          7.565658,\n"
		"          57.961921\n"
		"        ]\n"
		"      }\n"
		"    }\n"
		"  ]\n"
		"}";

	auto obj = createFromInput(json, nullptr);
	ASSERT_TRUE(obj != nullptr);

	auto jsonObj = nn_dynamic_pointer_cast<GeoJson>(obj);
	ASSERT_TRUE(jsonObj != nullptr);

	EXPECT_EQ(jsonObj->exportToJSON(&(JSONFormatter::create()->setSchema("foo"))), json);
}
