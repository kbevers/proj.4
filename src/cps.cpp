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

// TODO: Move to another class
struct COMMONPOINTS *cp_init(projCtx ctx, struct projFileAPI_t* fileapi)
{
	PAFile fid = (PAFile)fileapi;
	struct COMMONPOINTS *cp;
	int id_end;
	char header[160];

	cp = (struct COMMONPOINTS *) pj_malloc(sizeof(struct COMMONPOINTS));

	if (pj_ctx_fread(ctx, header, sizeof(header), 1, fid) != 1)
	{
		pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_CPL);
		return nullptr;
	}

	//cp ->cvs = nullptr; ?

	return cp;

}

PJ_COMMONPOINTS *pj_commonpoints_init(projCtx ctx, const char *cp_name)
{
	PJ_COMMONPOINTS *commonPoints;
	PAFile fp;
	size_t header_size = 0;

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

	// TODO: Sjekk header
	struct COMMONPOINTS *cp = cp_init( ctx, (struct projFileAPI_t*)fp );

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

	if (gi->cp->pJ_LP_PairList != nullptr)
	{
		pj_release_lock();
		return 1;
	}

	memcpy(&cp_tmp, gi->cp, sizeof(struct COMMONPOINTS));
	 
	//if (strcmp(gi->format, "ctable") == 0)

	PAFile fid;
	int result;
	
	fid = pj_open_lib(ctx, gi->filename, "rb");

	if (fid == nullptr)
	{
		pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_CPL);
		pj_release_lock();
		return 0;
	}

	size_t a_size;	

	if (pj_ctx_fread(ctx, &cp_tmp, sizeof(COMMONPOINTS), a_size, fid) != a_size)
	{

	};

	pj_ctx_fclose(ctx, fid);

	//if (gi->cp->pJ_LP_PairList != nullptr)

	return 1;
}
