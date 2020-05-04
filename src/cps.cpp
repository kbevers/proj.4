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
#include "cplist.hpp"

#include <algorithm>
#include <cmath>

NS_PROJ_START

LPZ_Pair::LPZ_Pair() = default;

// ---------------------------------------------------------------------------

Common_Points::Common_Points() = default;

// ---------------------------------------------------------------------------

Common_Points::Common_Points(std::unique_ptr<File> &&fp, const std::string &nameIn, const std::string &format, int noOfPoints) 
	: m_fp(std::move(fp)), m_name(nameIn), m_format(format), m_noOfPoints(noOfPoints)
{	
}

// ---------------------------------------------------------------------------

Common_Points::~Common_Points() = default;

// ---------------------------------------------------------------------------
 
Common_Points *Common_Points::open(PJ_CONTEXT *ctx, std::unique_ptr<File> fp, const std::string &filename)
{
	unsigned char header[160];

	if (fp->read(header, sizeof(header)) != sizeof(header))
	{
		pj_ctx_set_errno(ctx, PJD_ERR_FAILED_TO_LOAD_CPT);
		return nullptr;
	}

	__int32 noOfPoints;
	memcpy(&noOfPoints, header + 0, 4);	

	if (noOfPoints < 4)
		return nullptr;
	
	// TODO: Add name in CPT-file.
	// TODO: Add licence in CPT-file.
	 std::string name = filename;
	 std::string format = "cpt";
	//memcpy(&name, header + 4, 8);

	return new Common_Points(std::move(fp), name, format, noOfPoints);
}

// ---------------------------------------------------------------------------

bool Common_Points::load(PJ_CONTEXT *ctx)
{
	if (m_LpzPairList.size() == NoOfPoints())
		return true;

	auto pointPair = new LPZ_Pair();

	unsigned long offset = 4;
	m_fp->seek(offset);

	while (m_fp->read(pointPair, sizeof(LPZ_Pair) - offset) == sizeof(LPZ_Pair) - offset)
	{
		m_LpzPairList.push_back(*pointPair);
		pointPair = new LPZ_Pair();
	}

	if (m_LpzPairList.size() != NoOfPoints())
		return false;

	return true;
}

const Common_Points *Common_Points::cpAt(double lon, double lat) const
{
	double coslat = cos(lat);
	 
	for (auto&& pair : m_LpzPairList)
	{
		auto point = pair.FromPoint();
		double deltaPhi = point.phi - lat;
		double deltaLam = (point.lam - lon) * coslat;

		if (hypot(deltaPhi, deltaLam) < 1.0)
			return this;
	}	
	return nullptr;
}; 

NS_PROJ_END
