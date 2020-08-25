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
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <cmath>

#include "proj_internal.h"
#include "pps.hpp"
#include "pps_set.hpp"
#include "geojsonPolygon.hpp"

NS_PROJ_START

LPZ_Pair::LPZ_Pair() = default;

// ---------------------------------------------------------------------------

PointPairs::PointPairs() = default;

// ---------------------------------------------------------------------------

PointPairs::PointPairs(std::unique_ptr<File> &&fp, const std::string &nameIn, const std::string &format)
	: m_fp(std::move(fp)), m_name(nameIn), m_format(format)
{
}

// ---------------------------------------------------------------------------

PointPairs::~PointPairs() = default;

// ---------------------------------------------------------------------------
 
PointPairs *PointPairs::open(PJ_CONTEXT *ctx, std::unique_ptr<File> fp, const std::string &filename)
{
	unsigned char header[404];

	if (fp->read(header, sizeof(header)) != sizeof(header))
	{
		pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return nullptr;
	}

	int noOfPoints = 0;
	memcpy(&noOfPoints, header + 400, 4);	

	if (noOfPoints < 4)
		return nullptr;
	
	std::string name = filename;
	std::string format = "cpt";

	return new PointPairs(std::move(fp), name, format);
}

// ---------------------------------------------------------------------------

bool PointPairs::load(PJ_CONTEXT *ctx)
{
	if ((int)m_LpzPairList.size() > 0 && (int)m_LpzPairList.size() == NoOfPoints())
		return true;
	
	unsigned long offset = 404;
	unsigned long asize = 4;

	if (!m_fp->seek(offset))
		return false;

	auto pointPair = new LPZ_Pair();
	
	while (m_fp->read(pointPair = new LPZ_Pair(), sizeof(LPZ_Pair) - asize) == sizeof(LPZ_Pair) - asize)
	{
		m_LpzPairList.push_back(*pointPair);
		pointPair = new LPZ_Pair();
	}

	if ((int)m_LpzPairList.size() != NoOfPoints())
	{
		pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return false;
	}
	return true;
}

// THIS IS A TEST
bool PointPairs::loadGeoJson(PJ_CONTEXT *ctx)
{
	if (NoOfPoints() > 0)
		return true;
	 
	auto geoJsonSource = geoJson::GeoJson::openGeoJson(ctx, "EUREF89_NGO48_20081014_source.geojson");
	if (!geoJsonSource)
		return false;

	auto geoJsonTarget = geoJson::GeoJson::openGeoJson(ctx, "EUREF89_NGO48_20081014_target.geojson");
	if (!geoJsonTarget)
		return false;

	auto pointPairSets = new PointPairsSet();

	for (auto &feature : geoJsonSource->featuresMap())
	{
		auto name = feature.first;
		auto featureSource = feature.second;

		if (geoJsonTarget->FeatureExits(name))
		{
			auto featureTarget = geoJsonTarget->GetFeature(name);

			double xSource = featureSource->Point()->X_rad();
			double ySource = featureSource->Point()->Y_rad();

			double xTarget = featureTarget->Point()->X_rad();
			double yTarget = featureTarget->Point()->Y_rad();

			int areaId = featureTarget->AreaId();
			auto name = featureTarget->Name();
		
			auto pair = new LPZ_Pair();

			pair->SetFromPointPosition(xSource, ySource);
			pair->SetToPointPosition(xTarget, yTarget);
			pair->Area(areaId);
			pair->Name(strdup(name.c_str())); 

			m_LpzPairList.push_back(*pair);
		}
	}
	return true;
}

const PointPairs *PointPairs::pairsAt(double lon, double lat, double maxdist) const
{
	double coslat = cos(lat);

	for (auto&& pair : m_LpzPairList)
	{
		auto point = pair.FromPoint();

		double deltaPhi = point.phi - lat;
		double deltaLam = (point.lam - lon) * coslat;

		if (hypot(deltaPhi, deltaLam) < maxdist)
			return this;
	}
	return nullptr;
};

NS_PROJ_END
