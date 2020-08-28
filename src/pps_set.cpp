/************************************************************************
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
***********************************************************************/
#define PJ_LIB__

#ifndef FROM_PROJ_CPP
#define FROM_PROJ_CPP
#endif

#include <errno.h>
#include <stddef.h>
#include <string.h>

#include "pps.hpp"
#include "pps_set.hpp"
#include "geojsonPolygon.hpp"
#include "filemanager.hpp"
#include "proj/internal/internal.hpp"
#include "proj.h"
#include "proj_internal.h"

#include <algorithm>
#include <cmath>

#define PJ_MAX_PATH_LENGTH 1024

NS_PROJ_START

using namespace internal;

PointPairsSet::PointPairsSet() = default;

// ---------------------------------------------------------------------------

PointPairsSet::~PointPairsSet() = default;

// ---------------------------------------------------------------------------

std::unique_ptr<PointPairsSet> PointPairsSet::open(PJ_CONTEXT *ctx, const std::string &filename)
{
	if (filename == "null")
	{
		auto set = std::unique_ptr<PointPairsSet>(new PointPairsSet());
		set->m_name = filename;
		set->m_format = "null";
	 	set->m_pairs.push_back(std::unique_ptr<PointPairs>(new PointPairs()));

		return set;
	}
	auto fp = FileManager::open_resource_file(ctx, filename.c_str());
	if (!fp)
		return nullptr;
	
	const auto actualName(fp->name());

	if (ends_with(actualName, "cpt") || ends_with(actualName, "CPT"))
	{
		auto pp = PointPairs::open(ctx, std::move(fp), filename);
	
		if (!pp)
			return nullptr;
		
		if (!pp->load(ctx))
			return nullptr;
		
		auto set = std::unique_ptr<PointPairsSet>(new PointPairsSet());
		set->m_name = filename;
		set->m_format = "cpt";
		set->m_pairs.push_back(std::unique_ptr<PointPairs>(pp));

		return set;
	}
	return nullptr;
};

std::unique_ptr<PointPairsSet> PointPairsSet::open(PJ_CONTEXT *ctx, const std::string &sourcename, const std::string &targetname)
{
	if (sourcename == "null" || targetname == "null")
	{
		auto set = std::unique_ptr<PointPairsSet>(new PointPairsSet());

		set->m_sourceName = sourcename;
		set->m_targetName = targetname;
		set->m_format = "null";
		set->m_pairs.push_back(std::unique_ptr<PointPairs>(new PointPairs()));

		return set;
	}

	auto fpSource = FileManager::open_resource_file(ctx, sourcename.c_str());
	if (!fpSource)
		return nullptr;

	auto fpTarget = FileManager::open_resource_file(ctx, targetname.c_str());
	if (!fpTarget)
		return nullptr;

	auto geoJsonSource = geoJson::GeoJson::openGeoJson(ctx, sourcename);
	if (!geoJsonSource)
		return nullptr;

	auto geoJsonTarget = geoJson::GeoJson::openGeoJson(ctx, targetname);
	if (!geoJsonTarget)
		return nullptr;

	auto pointPairs = new PointPairs();
	 
	// Pairing source and target points
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
			auto nameStr = featureTarget->Name();

			auto pair = new LPZ_Pair();

			pair->Area(areaId);
			pair->Name(strdup(nameStr.c_str()));
			pair->SetFromPointPosition(xSource, ySource);
			pair->SetToPointPosition(xTarget, yTarget);
			 
			pointPairs->LpzPairList().push_back(*pair);
		}
	}
	
	auto set = std::unique_ptr<PointPairsSet>(new PointPairsSet());
	set->m_sourceName = sourcename;
	set->m_targetName = targetname;
	set->m_format = "geojson";
	set->m_pairs.push_back(std::unique_ptr<PointPairs>(pointPairs));
	
	return set;
}

// ---------------------------------------------------------------------------

ListOfPpSet pj_pp_init(PJ *P, const char *ppkey)
{
	std::string key("s");
	key += ppkey;

	const char *ppnames = pj_param(P->ctx, P->params, key.c_str()).s;
	if (ppnames == nullptr)
		return {};
	
	auto list = internal::split(std::string(ppnames), ',');
	ListOfPpSet pps;

	for (const auto &ppnameStr : list)
	{
		const char *ppname = ppnameStr.c_str();
		
		if (ppname[0] == '@')
			ppname++;
	
		auto ppSet = PointPairsSet::open(P->ctx, ppname);
		 
		if (!ppSet)
		{
			pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
			
			return {};
		}
		else
			pps.emplace_back(std::move(ppSet));
	}
	return pps;
} 

// ---------------------------------------------------------------------------

ListOfPpSet pj_pp_init(PJ *P, const char *sourcekey, const char *targetkey)
{
	std::string keyS("s");
	std::string keyT("s");

	keyS += sourcekey;
	keyT += targetkey;

	ListOfPpSet pps;

	const char *sourcename;
	const char *targetname;

	const char *sourcenames = pj_param(P->ctx, P->params, keyS.c_str()).s;
	if (sourcenames == nullptr)
		return {};
	
	const char *targetnames = pj_param(P->ctx, P->params, keyT.c_str()).s;
	if (targetnames == nullptr)
		return {};

	auto listSource = internal::split(std::string(sourcenames), ',');

	for (const auto &sourceStr : listSource)
		sourcename = sourceStr.c_str();

	auto listTarget = internal::split(std::string(targetnames), ',');
	
	for (const auto &targetStr : listTarget)
	    targetname = targetStr.c_str();

	auto ppSet = PointPairsSet::open(P->ctx, sourcename, targetname);
  
	if (!ppSet)
	{
		pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_GEOJSON);

		return {};
	}
	else
		pps.emplace_back(std::move(ppSet));

	return pps;
}

// ---------------------------------------------------------------------------

PointPairs* findPointPairs(PJ *P, const ListOfPpSet &pps, const PJ_LPZ &input, double maxdist)
{
	// Scaled km to radians:
	maxdist /= P->a / 1000.0;

	for (const auto &ppSet : pps)
	{
		if (ppSet->pairsAt(input.lam, input.phi, maxdist) != nullptr)
			return ppSet->pairsAt(input.lam, input.phi, maxdist);
	}
	return nullptr;
}

// ---------------------------------------------------------------------------

PointPairs *PointPairsSet::pairsAt(double lon, double lat, double maxdist) const
{
	for (const auto &pairs : m_pairs)
	{
		if (pairs->pairsAt(lon, lat, maxdist) != nullptr)
			return pairs.get();
	}
	return nullptr;
}

NS_PROJ_END
 