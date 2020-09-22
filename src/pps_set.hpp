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

//#ifndef CPS_HPP_INCLUDED
//#define CPS_HPP_INCLUDED

#include <memory>
#include <vector>

#include "pps.hpp"
#include "proj.h"
#include "proj/util.hpp"

NS_PROJ_START

// ---------------------------------------------------------------------------

class PROJ_GCC_DLL PointPairsSet
{
private:
protected:
	std::string m_name{};
	std::string m_sourceName{};
	std::string m_targetName{};
	std::string m_format{};
	std::vector<std::unique_ptr<PointPairs>> m_pairs {};
public:
	PROJ_FOR_TEST PointPairsSet();
	PROJ_FOR_TEST virtual ~PointPairsSet();
	PROJ_FOR_TEST static std::unique_ptr<PointPairsSet> open(PJ_CONTEXT *ctx, const std::string &filename);
	PROJ_FOR_TEST static std::unique_ptr<PointPairsSet> open(PJ_CONTEXT *ctx, const std::string &sourcename, const std::string &targetname);
	PROJ_FOR_TEST const std::string &Name() const { return m_name; }
	PROJ_FOR_TEST void Name(std::string name) { m_name = std::move(name); }
	PROJ_FOR_TEST const std::string &SourceName() const { return m_sourceName; }
	PROJ_FOR_TEST void SourceName(std::string name) { m_sourceName = std::move(name); }
	PROJ_FOR_TEST const std::string &TargetName() const { return m_targetName; }
	PROJ_FOR_TEST void TargetName(std::string name) { m_targetName = std::move(name); }
	PROJ_FOR_TEST const std::string &Format() const { return m_format; }
	PROJ_FOR_TEST void Format(std::string format) { m_format = std::move(format); }
	PROJ_FOR_TEST const std::vector<std::unique_ptr<PointPairs>> &Pairs() const { return m_pairs; }
	PROJ_FOR_TEST PointPairs *pairsAt(double lon, double lat, double maxdist) const;
};

// ---------------------------------------------------------------------------

typedef std::vector<std::unique_ptr<PointPairsSet>> ListOfPpSet;

// ---------------------------------------------------------------------------

ListOfPpSet pj_pp_init(PJ *P, const char *ppkey);

// ---------------------------------------------------------------------------

ListOfPpSet pj_pp_init(PJ *P, const char *sourcekey, const char *targetkey);

// ---------------------------------------------------------------------------

PointPairs *findPointPairs(PJ *P, const ListOfPpSet &cps, const PJ_LPZ &input, double maxdist = 0.1);

NS_PROJ_END

//#endif