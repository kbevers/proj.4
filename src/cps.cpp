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

#include "proj_internal.h"
#include "cps.hpp"

#include <algorithm>
#include <cmath>

NS_PROJ_START

Common_Points::Common_Points() = default;

// ---------------------------------------------------------------------------

Common_Points::~Common_Points() = default;

// ---------------------------------------------------------------------------

std::unique_ptr<Common_Points> *Common_Points::open(PJ_CONTEXT *ctx, std::unique_ptr<File> fp, const std::string &filename)
{
	unsigned char header[160];

	if (fp->read(header, sizeof(header)) != sizeof(header))
	{
		pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return nullptr;
	}


	return nullptr;
}

NS_PROJ_END

// TODO: Move to another class
struct COMMONPOINTS *cp_init(projCtx ctx, struct projFileAPI_t* fileapi)
{
	PAFile fid = (PAFile)fileapi;
	struct COMMONPOINTS *cp;	 
	char header[160];

	cp = (struct COMMONPOINTS *) pj_malloc(sizeof(struct COMMONPOINTS));

	if (pj_ctx_fread(ctx, header, sizeof(header), 1, fid) != 1)
	{
		pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return nullptr;
	}

	cp->noOfPoints = 0;
	cp->pJ_LPZ_PairList = nullptr;

	return cp;
}

PJ_COMMONPOINTS *pj_commonpoints_init(projCtx ctx, const char *cp_name)
{
	PJ_COMMONPOINTS *commonPoints;
	PAFile fp;	 

	errno = pj_errno = 0;
	ctx->last_errno = 0;

	commonPoints = (PJ_COMMONPOINTS *)pj_calloc(1, sizeof(PJ_COMMONPOINTS));

	if (!commonPoints)
	{
		pj_ctx_set_errno(ctx, ENOMEM);
		return nullptr;
	}

	commonPoints->cp_name = pj_strdup(cp_name);
	if (!commonPoints->cp_name)
	{
		pj_dalloc(commonPoints);
		pj_ctx_set_errno(ctx, ENOMEM);
		return nullptr;
	}
	commonPoints->filename = pj_strdup(cp_name);
	commonPoints->format = "missing";
	commonPoints->cp = nullptr;
	commonPoints->next = nullptr;

	if (!(fp = pj_open_lib(ctx, cp_name, "rb")))
	{
		ctx->last_errno = 0;
		return commonPoints;
	}

	if (!commonPoints->filename)
	{
		pj_dalloc(commonPoints->cp_name);
		pj_dalloc(commonPoints);
		pj_ctx_set_errno(ctx, ENOMEM);
		return nullptr;
	}
	 
	struct COMMONPOINTS *cp = cp_init(ctx, (struct projFileAPI_t*)fp);

	commonPoints->format = "commonpoints";
	commonPoints->cp = cp;

	if (cp == nullptr)
	{
		pj_log(ctx, PJ_LOG_DEBUG_MAJOR, "COMMONPOINTS cp is NULL.");
	}
	else
	{
		// TODO: Logging
		/*
		  pj_log( ctx, PJ_LOG_DEBUG_MAJOR,
                    "Ctable2 %s %dx%d: LL=(%.9g,%.9g) UR=(%.9g,%.9g)",
                    ct->id,
                    ct->lim.lam, ct->lim.phi,
                    ct->ll.lam * RAD_TO_DEG, ct->ll.phi * RAD_TO_DEG,
                    (ct->ll.lam + (ct->lim.lam-1)*ct->del.lam) * RAD_TO_DEG,
                    (ct->ll.phi + (ct->lim.phi-1)*ct->del.phi) * RAD_TO_DEG );
		
		*/
	}

	pj_ctx_fclose(ctx, fp);

	return commonPoints;
}

int pj_cp_load(projCtx_t* ctx, PJ_COMMONPOINTS *gi)
{
	struct COMMONPOINTS cp_tmp;

	if (gi == nullptr || gi->cp == nullptr)
		return 0;

	pj_acquire_lock();

	if (gi->cp->pJ_LPZ_PairList != nullptr || gi->cp->noOfPoints > 0)
	{
		pj_release_lock();
		return 1;
	}

	memcpy(&cp_tmp, gi->cp, sizeof(struct COMMONPOINTS));

	PAFile fid;
 	
	fid = pj_open_lib(ctx, gi->filename, "rb");

	if (fid == nullptr)
	{
		pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		pj_release_lock();

		return 0;
	}

	std::vector<PJ_LPZ_Pair> *pJLPZList = cp_tmp.pJ_LPZ_PairList;
	pJLPZList = (std::vector<PJ_LPZ_Pair> *)pj_calloc(1, sizeof(struct std::vector<PJ_LPZ_Pair>));

	int noOfPoints = 0;
	long a_size = sizeof(struct PJ_LPZ_Pair) - 4;

	pj_ctx_fread(ctx, &noOfPoints, sizeof(__int32), 1, fid);

	for (int i = 0; i < noOfPoints; i++)
	{
		pj_ctx_fseek(ctx, fid, a_size * i + sizeof(__int32), SEEK_SET);
		PJ_LPZ_Pair *pair = (PJ_LPZ_Pair *)pj_calloc(1, sizeof(struct PJ_LPZ_Pair));
		pj_ctx_fread(ctx, pair, sizeof(struct PJ_LPZ_Pair), 1, fid);

		pJLPZList->push_back(*pair);
	}

	gi->cp->noOfPoints = noOfPoints;
	gi->cp->pJ_LPZ_PairList = pJLPZList;

	pj_ctx_fclose(ctx, fid);

	return 1;
}
