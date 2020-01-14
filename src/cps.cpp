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

	if (!(fp = pj_open_lib(ctx, cp_name, "rb")))
	{
		ctx->last_errno = 0;
		return commonPoints;
	}
	commonPoints->filename = pj_strdup(cp_name);

	if (!commonPoints->filename)
	{
		pj_dalloc(commonPoints->cp_name);
		pj_dalloc(commonPoints);
		pj_ctx_set_errno(ctx, ENOMEM);

		return nullptr;
	}

	// TODO: Impl. here.
	// Ala:   struct CTABLE *ct = nad_ctable2_init( ctx, (struct projFileAPI_t*)fp );



	pj_ctx_fclose(ctx, fp);

	return commonPoints;
}

// TODO: Move to another class
struct COMMONPOINTS *cp_init(projCtx ctx, struct projFileAPI_t* fileapi)
{
	PAFile fid = (PAFile)fileapi;
	struct COMMONPOINTS *cp;
	int id_end;


	cp = (struct COMMONPOINTS *) pj_malloc(sizeof(struct COMMONPOINTS));

	return cp;

}
