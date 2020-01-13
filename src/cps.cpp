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

PJ_COMMONPOINTS *pj_commonpoints_init(projCtx ctx, const char *fileName)
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

	commonPoints->filename = pj_strdup(fileName);
	if (!commonPoints->filename)
	{
		pj_dalloc(commonPoints);
		pj_ctx_set_errno(ctx, ENOMEM);
		return nullptr;
	}

	if (!(fp = pj_open_lib(ctx, fileName, "rb")))
	{
		ctx->last_errno = 0; /* don't treat as a persistent error */
		return commonPoints;
	}
	commonPoints->filename = pj_strdup(fileName);

	if (!commonPoints->filename)
	{
		pj_dalloc(commonPoints->filename);
		pj_dalloc(commonPoints);
		pj_ctx_set_errno(ctx, ENOMEM);
		return nullptr;
	}

	pj_ctx_fclose(ctx, fp);

	return commonPoints;
}

