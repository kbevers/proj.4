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

// ---------------------------------------------------------------------------

ListOfPpSet pj_cp_init(PJ *P, const char *cpkey)
{
	std::string key("s");
	key += cpkey;

	const char *ppnames = pj_param(P->ctx, P->params, key.c_str()).s;
	if (ppnames == nullptr)
		return {};
	
	auto list = internal::split(std::string(ppnames), ',');
	ListOfPpSet cps;

	for (const auto &ppnameStr : list)
	{
		const char *ppname = ppnameStr.c_str();
		bool canFail = false;
		if (ppname[0] == '@')
		{
			canFail = true;
			ppname++;
		}
		auto ppSet = PointPairsSet::open(P->ctx, ppname);
		if (!ppSet)
		{
			if (!canFail)
			{
		 	//	if (proj_context_errno(P->ctx) != PJD_ERR_NETWORK_ERROR) 
				{
					pj_ctx_set_errno(P->ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
				}
				return {};
			}
			pj_ctx_set_errno(P->ctx, 0);
		}
		else
		{
			cps.emplace_back(std::move(ppSet));
		}	
	}
	return cps;
} 

// ---------------------------------------------------------------------------

PointPairs* findPointPairs(const ListOfPpSet &pps, const PJ_LPZ &input)
{
	for (const auto &ppSet : pps)
	{
		if (ppSet->pairsAt(input.phi, input.phi) != nullptr)
		{
			return ppSet->pairsAt(input.phi, input.phi);
		}
		//return nullptr;
		//	return cpSet;

		/*
		cpSetOut = cpSet.get();

		if (cpSetOut == nullptr)
			return nullptr;

		if (cpSetOut->Cps().size() == 0)
			return nullptr;
		*/
		// TODO: Add extent area in cpt-file.
		//	return cpSet->Cps();
	}
	return nullptr;
}

// ---------------------------------------------------------------------------

PointPairs *PointPairsSet::pairsAt(double lon, double lat) const
{
	for (const auto &pairs : m_pairs)
	{ 
		if (pairs->ppAt(lon, lat) != nullptr)
			return pairs.get();
	}
	return nullptr;
}

NS_PROJ_END
 