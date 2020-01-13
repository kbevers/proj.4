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
#include <stddef.h>
#include <string.h>

#include "proj.h"
#include "proj_internal.h"

static PJ_COMMONPOINTS *list = nullptr;
#define PJ_MAX_PATH_LENGTH 1024

static int pj_cplist_merge(projCtx ctx, const char *cp_name, PJ_COMMONPOINTS ***p_list, int *p_listcount, int *p_list_max)
{
	int match = 0;
	PJ_COMMONPOINTS *this_cp, *tail = nullptr;
	
	for (this_cp = list; this_cp != nullptr; this_cp = this_cp->next)
	{
		if (strcmp(this_cp->cp_name, cp_name) == 0)
		{
			match = 1;

			if (this_cp->cp == nullptr)
				return 0;

			if (*p_listcount >= *p_list_max - 2)
			{
				PJ_COMMONPOINTS **new_list;

				int new_max = *p_list_max + 20;

				new_list = (PJ_COMMONPOINTS **)pj_calloc(new_max, sizeof(void *));
				
				if (!new_list) 
				{
					pj_ctx_set_errno(ctx, ENOMEM);
					return 0;
				}
				if (*p_list != nullptr)
				{
					memcpy(new_list, *p_list, sizeof(void *) * (*p_list_max));
					pj_dalloc(*p_list);
				}

			}


		}
	}


	// TODO: Recursive...

	return 0;
}

PJ_COMMONPOINTS **pj_cplist(projCtx ctx, const char *lists, int *list_count)
{
	const char *s;
	PJ_COMMONPOINTS **list = nullptr;
	int list_max = 0;

	pj_errno = 0;
	*list_count = 0;
	 
	for (s = lists; *s != '\0'; )
	{
		size_t end_char;
		int required = 1;
		char   name[PJ_MAX_PATH_LENGTH];

		if (*s == '@')
		{
			required = 0;
			s++;
		}

		for (end_char = 0; s[end_char] != '\0' && s[end_char] != ','; end_char++)
		{
		}

		if (end_char >= sizeof(name))
		{
			pj_dalloc(list);
			pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_CPL);
			pj_release_lock();
			return nullptr;
		}

		strncpy(name, s, end_char);
		name[end_char] = '\0';

		s += end_char;
		if (*s == ',')
			s++;

		if (!pj_cplist_merge(ctx, name, &list, list_count, &list_max) && required)
		{
			pj_dalloc(list);
			pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_CPL);
			pj_release_lock();
			return nullptr;
		}
		else
			pj_errno = 0;
	}
	pj_release_lock();

	return list;
}
