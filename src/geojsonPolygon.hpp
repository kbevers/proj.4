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

#include <memory>
#include <vector>

//#include "proj.h"
//#include "proj_internal.h"
//#include "proj/util.hpp"
#include "point_in_polygon.h"

namespace
{
	struct geoJsonMultiPolygon
	{
		std::string name;
		vector<PolygonPoint> *pointList;
	};
}

typedef std::vector<std::unique_ptr<geoJsonMultiPolygon>> ListOfMultiPolygon;

ListOfMultiPolygon pj_polygon_init(PJ *P, const char *polygons);
 
void testReadGeojson(/*char* fileName*/);
bool pointIsInArea(PJ_LP pointPJ_LP, char* fileName);
int areaIdPoint(PJ_LP *lp);