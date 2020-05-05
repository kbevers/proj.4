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

#ifndef CPS_HPP_INCLUDED
#define CPS_HPP_INCLUDED

#include <memory>
#include <vector>

#include "proj.h"
#include "proj/util.hpp"
#include "filemanager.hpp"

NS_PROJ_START

class PROJ_GCC_DLL LPZ_Pair
{
private:
	char name[8];
	PJ_LPZ m_fromPoint;
	PJ_LPZ m_toPoint;
	__int32 m_area = 0;
	double m_dist = 0.0;
protected:
public:
	PROJ_FOR_TEST LPZ_Pair();
	PROJ_FOR_TEST const PJ_LPZ &FromPoint() const { return m_fromPoint; }
	PROJ_FOR_TEST const PJ_LPZ &ToPoint() const { return m_toPoint; }
	PROJ_FOR_TEST const __int32 &Area() const { return m_area; }
	PROJ_FOR_TEST const double Distance() const { return m_dist; }
	PROJ_FOR_TEST void SetDistance(double dist) { m_dist = dist; }
};

// ---------------------------------------------------------------------------

class PROJ_GCC_DLL PointPairs
{
private:
protected:
	int m_noOfPoints = 0;
	std::string m_name { };
	std::string m_format { };
	std::vector<LPZ_Pair> m_LpzPairList { };
	
	PJ_CONTEXT *m_ctx;
	std::unique_ptr<File> m_fp;
public:
	PROJ_FOR_TEST PointPairs();
	PROJ_FOR_TEST PointPairs(std::unique_ptr<File> &&fp, const std::string &nameIn, const std::string &format, int noOfPoints);
    PROJ_FOR_TEST virtual ~PointPairs();
	PROJ_FOR_TEST int NoOfPoints() const { return m_noOfPoints; }
	PROJ_FOR_TEST const std::string &Name() const { return m_name; }
	PROJ_FOR_TEST const std::string &Format() const { return m_format; }
	PROJ_FOR_TEST const std::vector<LPZ_Pair> &LpzPairList() const { return m_LpzPairList; }
    PROJ_FOR_TEST static PointPairs *open(PJ_CONTEXT *ctx, std::unique_ptr<File> fp, const std::string &filename);
	PROJ_FOR_TEST bool load(PJ_CONTEXT *ctx);
	PROJ_FOR_TEST const PointPairs *pairsAt(double lon, double lat, double mindist = 0.1) const;
};
NS_PROJ_END

#endif