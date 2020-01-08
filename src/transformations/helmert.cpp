/***********************************************************************

             3-, 4-and 7-parameter shifts, and their 6-, 8-
                and 14-parameter kinematic counterparts.

                    Thomas Knudsen, 2016-05-24

************************************************************************

  Implements 3(6)-, 4(8) and 7(14)-parameter Helmert transformations for
  3D data.

  Also incorporates Molodensky-Badekas variant of 7-parameter Helmert
  transformation, where the rotation is not applied regarding the centre
  of the spheroid, but given a reference point.

  Primarily useful for implementation of datum shifts in transformation
  pipelines.

************************************************************************

Thomas Knudsen, thokn@sdfe.dk, 2016-05-24/06-05
Kristian Evers, kreve@sdfe.dk, 2017-05-01
Even Rouault, even.roault@spatialys.com
Last update: 2018-10-26

************************************************************************
* Copyright (c) 2016, Thomas Knudsen / SDFE
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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstddef>
#include <algorithm>

#include "proj_internal.h"
#include "geocent.h"
#include "point_in_polygon.h"
#include "proj\internal\nlohmann\json.hpp"
#include "Eigen\Eigen"

using namespace Eigen; 
using json = nlohmann::json;

PROJ_HEAD(helmert, "3(6)-, 4(8)- and 7(14)-parameter Helmert shift");
PROJ_HEAD(molobadekas, "Molodensky-Badekas transformation");

static PJ_XYZ helmert_forward_3d (PJ_LPZ lpz, PJ *P);
static PJ_LPZ helmert_reverse_3d (PJ_XYZ xyz, PJ *P);

/***********************************************************************/
namespace { // anonymous namespace
struct pj_opaque_helmert {
/************************************************************************
    Projection specific elements for the "helmert" PJ object
************************************************************************/
    PJ_XYZ xyz;
    PJ_XYZ xyz_0;
    PJ_XYZ dxyz;
    PJ_XYZ refp;
    PJ_OPK opk;
    PJ_OPK opk_0;
    PJ_OPK dopk;
    double scale;
    double scale_0;
    double dscale;
    double theta;
    double theta_0;
    double dtheta;
    double R[3][3];
    double t_epoch, t_obs;
    int no_rotation, exact, fourparam;
    int is_position_vector; /* 1 = position_vector, 0 = coordinate_frame */
};
} // anonymous namespace

/* Make the maths of the rotation operations somewhat more readable and textbook like */
#define R00 (Q->R[0][0])
#define R01 (Q->R[0][1])
#define R02 (Q->R[0][2])

#define R10 (Q->R[1][0])
#define R11 (Q->R[1][1])
#define R12 (Q->R[1][2])

#define R20 (Q->R[2][0])
#define R21 (Q->R[2][1])
#define R22 (Q->R[2][2])

struct CommonPointPair
{
	std::string name;
	PJ_LP fromPoint;
	PJ_LP toPoint;
	__int32 area;
	double dist;
};
	
/**************************************************************************/
static void update_parameters(PJ *P) {
/***************************************************************************

    Update transformation parameters.
    ---------------------------------

    The 14-parameter Helmert transformation is at it's core the same as the
    7-parameter transformation, since the transformation parameters are
    projected forward or backwards in time via the rate of changes of the
    parameters. The transformation parameters are calculated for a specific
    epoch before the actual Helmert transformation is carried out.

    The transformation parameters are updated with the following equation [0]:

                      .
    P(t) = P(EPOCH) + P * (t - EPOCH)

                                                              .
    where EPOCH is the epoch indicated in the above table and P is the rate
    of that parameter.

    [0] http://itrf.ign.fr/doc_ITRF/Transfo-ITRF2008_ITRFs.txt

*******************************************************************************/

    struct pj_opaque_helmert *Q = (struct pj_opaque_helmert *) P->opaque;
    double dt = Q->t_obs - Q->t_epoch;

    Q->xyz.x = Q->xyz_0.x + Q->dxyz.x * dt;
    Q->xyz.y = Q->xyz_0.y + Q->dxyz.y * dt;
    Q->xyz.z = Q->xyz_0.z + Q->dxyz.z * dt;

    Q->opk.o = Q->opk_0.o + Q->dopk.o * dt;
    Q->opk.p = Q->opk_0.p + Q->dopk.p * dt;
    Q->opk.k = Q->opk_0.k + Q->dopk.k * dt;

    Q->scale = Q->scale_0 + Q->dscale * dt;

    Q->theta = Q->theta_0 + Q->dtheta * dt;

    /* debugging output */
    if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE) {
        proj_log_trace(P, "Transformation parameters for observation "
                       "t_obs=%g (t_epoch=%g):", Q->t_obs, Q->t_epoch);
        proj_log_trace(P, "x: %g", Q->xyz.x);
        proj_log_trace(P, "y: %g", Q->xyz.y);
        proj_log_trace(P, "z: %g", Q->xyz.z);
        proj_log_trace(P, "s: %g", Q->scale*1e-6);
        proj_log_trace(P, "rx: %g", Q->opk.o);
        proj_log_trace(P, "ry: %g", Q->opk.p);
        proj_log_trace(P, "rz: %g", Q->opk.k);
        proj_log_trace(P, "theta: %g", Q->theta);
    }
}

/**************************************************************************/
static void build_rot_matrix(PJ *P) {
/***************************************************************************

    Build rotation matrix.
    ----------------------

    Here we rename rotation indices from omega, phi, kappa (opk), to
    fi (i.e. phi), theta, psi (ftp), in order to reduce the mental agility
    needed to implement the expression for the rotation matrix derived over
    at https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions
    The relevant section is Euler angles ( z-’-x" intrinsic) -> Rotation matrix

    By default small angle approximations are used:
    The matrix elements are approximated by expanding the trigonometric
    functions to linear order (i.e. cos(x) = 1, sin(x) = x), and discarding
    products of second order.

    This was a useful hack when calculating by hand was the only option,
    but in general, today, should be avoided because:

    1. It does not save much computation time, as the rotation matrix
       is built only once and probably used many times (except when
       transforming spatio-temporal coordinates).

    2. The error induced may be too large for ultra high accuracy
       applications: the Earth is huge and the linear error is
       approximately the angular error multiplied by the Earth radius.

    However, in many cases the approximation is necessary, since it has
    been used historically: Rotation angles from older published datum
    shifts may actually be a least squares fit to the linearized rotation
    approximation, hence not being strictly valid for deriving the exact
    rotation matrix. In fact, most publicly available transformation
    parameters are based on the approximate Helmert transform, which is why
    we use that as the default setting, even though it is more correct to
    use the exact form of the equations.

    So in order to fit historically derived coordinates, the access to
    the approximate rotation matrix is necessary - at least in principle.

    Also, when using any published datum transformation information, one
    should always check which convention (exact or approximate rotation
    matrix) is expected, and whether the induced error for selecting
    the opposite convention is acceptable (which it often is).


    Sign conventions
    ----------------

    Take care: Two different sign conventions exist for the rotation terms.

    Conceptually they relate to whether we rotate the coordinate system
    or the "position vector" (the vector going from the coordinate system
    origin to the point being transformed, i.e. the point coordinates
    interpreted as vector coordinates).

    Switching between the "position vector" and "coordinate system"
    conventions is simply a matter of switching the sign of the rotation
    angles, which algebraically also translates into a transposition of
    the rotation matrix.

    Hence, as geodetic constants should preferably be referred to exactly
    as published, the "convention" option provides the ability to switch
    between the conventions.

***************************************************************************/
    struct pj_opaque_helmert *Q = (struct pj_opaque_helmert *) P->opaque;

    double  f,  t,  p;    /* phi/fi , theta, psi  */
    double cf, ct, cp;    /* cos (fi, theta, psi) */
    double sf, st, sp;    /* sin (fi, theta, psi) */

    /* rename   (omega, phi, kappa)   to   (fi, theta, psi)   */
    f = Q->opk.o;
    t = Q->opk.p;
    p = Q->opk.k;

    /* Those equations are given assuming coordinate frame convention. */
    /* For the position vector convention, we transpose the matrix just after. */
    if (Q->exact) {
        cf = cos(f);
        sf = sin(f);
        ct = cos(t);
        st = sin(t);
        cp = cos(p);
        sp = sin(p);


        R00 = ct*cp;
        R01 = cf*sp + sf*st*cp;
        R02 = sf*sp - cf*st*cp;

        R10 = -ct*sp;
        R11 =  cf*cp - sf*st*sp;
        R12 =  sf*cp + cf*st*sp;

        R20 =  st;
        R21 = -sf*ct;
        R22 =  cf*ct;
    } else{
        R00 =  1;
        R01 =  p;
        R02 = -t;

        R10 = -p;
        R11 =  1;
        R12 =  f;

        R20 =  t;
        R21 = -f;
        R22 =  1;
    }


    /*
        For comparison: Description from Engsager/Poder implementation
        in set_dtm_1.c (trlib)

        DATUM SHIFT:
        TO = scale * ROTZ * ROTY * ROTX * FROM + TRANSLA

             ( cz sz 0)         (cy 0 -sy)         (1   0  0)
        ROTZ=(-sz cz 0),   ROTY=(0  1   0),   ROTX=(0  cx sx)
             (  0  0 1)         (sy 0  cy)         (0 -sx cx)

        trp->r11  =  cos_ry*cos_rz;
        trp->r12  =  cos_rx*sin_rz + sin_rx*sin_ry*cos_rz;
        trp->r13  =  sin_rx*sin_rz - cos_rx*sin_ry*cos_rz;

        trp->r21  = -cos_ry*sin_rz;
        trp->r22  =  cos_rx*cos_rz - sin_rx*sin_ry*sin_rz;
        trp->r23  =  sin_rx*cos_rz + cos_rx*sin_ry*sin_rz;

        trp->r31  =  sin_ry;
        trp->r32  = -sin_rx*cos_ry;
        trp->r33  =  cos_rx*cos_ry;

        trp->scale = 1.0 + scale;
    */


    if (Q->is_position_vector) {
        double r;
        r = R01;    R01 = R10;    R10 = r;
        r = R02;    R02 = R20;    R20 = r;
        r = R12;    R12 = R21;    R21 = r;
    }

    /* some debugging output */
    if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_TRACE) {
        proj_log_trace(P, "Rotation Matrix:");
        proj_log_trace(P, "  | % 6.6g  % 6.6g  % 6.6g |", R00, R01, R02);
        proj_log_trace(P, "  | % 6.6g  % 6.6g  % 6.6g |", R10, R11, R12);
        proj_log_trace(P, "  | % 6.6g  % 6.6g  % 6.6g |", R20, R21, R22);
    }
}

/***********************************************************************/
static PJ_XY helmert_forward (PJ_LP lp, PJ *P) {
/***********************************************************************/
    struct pj_opaque_helmert *Q = (struct pj_opaque_helmert *) P->opaque;
    PJ_COORD point = {{0,0,0,0}};
    double x, y, cr, sr;
    point.lp = lp;

    cr = cos(Q->theta) * Q->scale;
    sr = sin(Q->theta) * Q->scale;
    x = point.xy.x;
    y = point.xy.y;

    point.xy.x =  cr*x + sr*y + Q->xyz_0.x;
    point.xy.y = -sr*x + cr*y + Q->xyz_0.y;

    return point.xy;
}

/***********************************************************************/
static PJ_LP helmert_reverse (PJ_XY xy, PJ *P) {
/***********************************************************************/
    struct pj_opaque_helmert *Q = (struct pj_opaque_helmert *) P->opaque;
    PJ_COORD point = {{0,0,0,0}};
    double x, y, sr, cr;
    point.xy = xy;

    cr = cos(Q->theta) / Q->scale;
    sr = sin(Q->theta) / Q->scale;
    x = point.xy.x - Q->xyz_0.x;
    y = point.xy.y - Q->xyz_0.y;

    point.xy.x =  x*cr - y*sr;
    point.xy.y =  x*sr + y*cr;

    return point.lp;
}

/***********************************************************************/
static PJ_XYZ helmert_forward_3d (PJ_LPZ lpz, PJ *P) {
/***********************************************************************/
    struct pj_opaque_helmert *Q = (struct pj_opaque_helmert *) P->opaque;
    PJ_COORD point = {{0,0,0,0}};
    double X, Y, Z, scale;

    point.lpz = lpz;

    if (Q->fourparam) {
        point.xy = helmert_forward(point.lp, P);
        return point.xyz;
    }

    if (Q->no_rotation) {
        point.xyz.x = lpz.lam + Q->xyz.x;
        point.xyz.y = lpz.phi + Q->xyz.y;
        point.xyz.z = lpz.z   + Q->xyz.z;
        return point.xyz;
    }

    scale = 1 + Q->scale * 1e-6;

    X = lpz.lam - Q->refp.x;
    Y = lpz.phi - Q->refp.y;
    Z = lpz.z - Q->refp.z;


    point.xyz.x = scale * ( R00 * X  +   R01 * Y   +   R02 * Z);
    point.xyz.y = scale * ( R10 * X  +   R11 * Y   +   R12 * Z);
    point.xyz.z = scale * ( R20 * X  +   R21 * Y   +   R22 * Z);

    point.xyz.x += Q->xyz.x; /* for Molodensky-Badekas, Q->xyz already incorporates the Q->refp offset */
    point.xyz.y += Q->xyz.y;
    point.xyz.z += Q->xyz.z;

    return point.xyz;
}

/***********************************************************************/
static PJ_LPZ helmert_reverse_3d (PJ_XYZ xyz, PJ *P) {
/***********************************************************************/
    struct pj_opaque_helmert *Q = (struct pj_opaque_helmert *) P->opaque;
    PJ_COORD point = {{0,0,0,0}};
    double X, Y, Z, scale;

    point.xyz = xyz;

    if (Q->fourparam) {
        point.lp = helmert_reverse(point.xy, P);
        return point.lpz;
    }

    if (Q->no_rotation) {
        point.xyz.x  =  xyz.x - Q->xyz.x;
        point.xyz.y  =  xyz.y - Q->xyz.y;
        point.xyz.z  =  xyz.z - Q->xyz.z;
        return point.lpz;
    }

    scale = 1 + Q->scale * 1e-6;

    /* Unscale and deoffset */
    X = (xyz.x - Q->xyz.x) / scale;
    Y = (xyz.y - Q->xyz.y) / scale;
    Z = (xyz.z - Q->xyz.z) / scale;

    /* Inverse rotation through transpose multiplication */
    point.xyz.x  =  ( R00 * X   +   R10 * Y   +   R20 * Z) + Q->refp.x;
    point.xyz.y  =  ( R01 * X   +   R11 * Y   +   R21 * Z) + Q->refp.y;
    point.xyz.z  =  ( R02 * X   +   R12 * Y   +   R22 * Z) + Q->refp.z;

    return point.lpz;
}

static PJ_COORD helmert_forward_4d (PJ_COORD point, PJ *P) {
    struct pj_opaque_helmert *Q = (struct pj_opaque_helmert *) P->opaque;

    /* We only need to rebuild the rotation matrix if the
     * observation time is different from the last call */
    double t_obs = (point.xyzt.t == HUGE_VAL) ? Q->t_epoch : point.xyzt.t;
    if (t_obs != Q->t_obs) {
        Q->t_obs = t_obs;
        update_parameters(P);
        build_rot_matrix(P);
    }

    point.xyz = helmert_forward_3d (point.lpz, P);

    return point;
}

static PJ_COORD helmert_reverse_4d (PJ_COORD point, PJ *P) {
    struct pj_opaque_helmert *Q = (struct pj_opaque_helmert *) P->opaque;

    /* We only need to rebuild the rotation matrix if the
     * observation time is different from the last call */
    double t_obs = (point.xyzt.t == HUGE_VAL) ? Q->t_epoch : point.xyzt.t;
    if (t_obs != Q->t_obs) {
        Q->t_obs = t_obs;
        update_parameters(P);
        build_rot_matrix(P);
    }

    point.lpz = helmert_reverse_3d (point.xyz, P);

    return point;
}

/* Arcsecond to radians */
#define ARCSEC_TO_RAD (DEG_TO_RAD / 3600.0)

static PJ* init_helmert_six_parameters(PJ* P) {
    struct pj_opaque_helmert *Q = static_cast<struct pj_opaque_helmert*>(pj_calloc (1, sizeof (struct pj_opaque_helmert)));
    if (nullptr==Q)
        return pj_default_destructor (P, ENOMEM);
    P->opaque = (void *) Q;

    /* In most cases, we work on 3D cartesian coordinates */
    P->left  = PJ_IO_UNITS_CARTESIAN;
    P->right = PJ_IO_UNITS_CARTESIAN;

    /* Translations */
    if (pj_param (P->ctx, P->params, "tx").i)
        Q->xyz_0.x = pj_param (P->ctx, P->params, "dx").f;

    if (pj_param (P->ctx, P->params, "ty").i)
        Q->xyz_0.y = pj_param (P->ctx, P->params, "dy").f;

    if (pj_param (P->ctx, P->params, "tz").i)
        Q->xyz_0.z = pj_param (P->ctx, P->params, "dz").f;

    /* Rotations */
    if (pj_param (P->ctx, P->params, "trx").i)
        Q->opk_0.o = pj_param (P->ctx, P->params, "drx").f * ARCSEC_TO_RAD;

    if (pj_param (P->ctx, P->params, "try").i)
        Q->opk_0.p = pj_param (P->ctx, P->params, "dry").f * ARCSEC_TO_RAD;

    if (pj_param (P->ctx, P->params, "trz").i)
        Q->opk_0.k = pj_param (P->ctx, P->params, "drz").f * ARCSEC_TO_RAD;

    /* Use small angle approximations? */
    if (pj_param (P->ctx, P->params, "bexact").i)
        Q->exact = 1;

    return P;
}

static PJ* read_convention(PJ* P) {

    struct pj_opaque_helmert *Q = (struct pj_opaque_helmert *)P->opaque;

    /* In case there are rotational terms, we require an explicit convention
     * to be provided. */
    if (!Q->no_rotation) {
        const char* convention = pj_param (P->ctx, P->params, "sconvention").s;
        if( !convention ) {
            proj_log_error (P, "helmert: missing 'convention' argument");
            return pj_default_destructor (P, PJD_ERR_MISSING_ARGS);
        }
        if( strcmp(convention, "position_vector") == 0 ) {
            Q->is_position_vector = 1;
        }
        else if( strcmp(convention, "coordinate_frame") == 0 ) {
            Q->is_position_vector = 0;
        }
        else {
            proj_log_error (P, "helmert: invalid value for 'convention' argument");
            return pj_default_destructor (P, PJD_ERR_INVALID_ARG);
        }

        /* historically towgs84 in PROJ has always been using position_vector
         * convention. Accepting coordinate_frame would be confusing. */
        if (pj_param_exists (P->params, "towgs84")) {
            if( !Q->is_position_vector ) {
                proj_log_error (P, "helmert: towgs84 should only be used with "
                                "convention=position_vector");
                return pj_default_destructor (P, PJD_ERR_INVALID_ARG);
            }
        }
    }

    return P;
}

/***********************************************************************/
PJ *TRANSFORMATION(helmert, 0) {
/***********************************************************************/

    struct pj_opaque_helmert *Q;

    if( !init_helmert_six_parameters(P) ) {
        return nullptr;
    }

    /* In the 2D case, the coordinates are projected */
    if (pj_param_exists (P->params, "theta")) {
        P->left  = PJ_IO_UNITS_PROJECTED;
        P->right = PJ_IO_UNITS_PROJECTED;
        P->fwd    = helmert_forward;
        P->inv    = helmert_reverse;
    }

    P->fwd4d  = helmert_forward_4d;
    P->inv4d  = helmert_reverse_4d;
    P->fwd3d  = helmert_forward_3d;
    P->inv3d  = helmert_reverse_3d;

    Q = (struct pj_opaque_helmert *)P->opaque;

    /* Detect obsolete transpose flag and error out if found */
    if (pj_param (P->ctx, P->params, "ttranspose").i) {
        proj_log_error (P, "helmert: 'transpose' argument is no longer valid. "
                        "Use convention=position_vector/coordinate_frame");
        return pj_default_destructor (P, PJD_ERR_INVALID_ARG);
    }

    /* Support the classic PROJ towgs84 parameter, but allow later overrides.*/
    /* Note that if towgs84 is specified, the datum_params array is set up   */
    /* for us automagically by the pj_datum_set call in pj_init_ctx */
    if (pj_param_exists (P->params, "towgs84")) {
        Q->xyz_0.x = P->datum_params[0];
        Q->xyz_0.y = P->datum_params[1];
        Q->xyz_0.z = P->datum_params[2];

        Q->opk_0.o = P->datum_params[3];
        Q->opk_0.p = P->datum_params[4];
        Q->opk_0.k = P->datum_params[5];

        /* We must undo conversion to absolute scale from pj_datum_set */
        if (0==P->datum_params[6])
            Q->scale_0 = 0;
        else
            Q->scale_0 = (P->datum_params[6] - 1) * 1e6;
    }

    if (pj_param (P->ctx, P->params, "ttheta").i) {
        Q->theta_0 = pj_param (P->ctx, P->params, "dtheta").f * ARCSEC_TO_RAD;
        Q->fourparam = 1;
        Q->scale_0 = 1.0; /* default scale for the 4-param shift */
    }

    /* Scale */
    if (pj_param (P->ctx, P->params, "ts").i) {
        Q->scale_0 = pj_param (P->ctx, P->params, "ds").f;
        if( Q->scale_0 <= -1.0e6 )
            return pj_default_destructor (P, PJD_ERR_INVALID_SCALE);
        if (pj_param (P->ctx, P->params, "ttheta").i && Q->scale_0 == 0.0)
            return pj_default_destructor (P, PJD_ERR_INVALID_SCALE);
    }

    /* Translation rates */
    if (pj_param(P->ctx, P->params, "tdx").i)
        Q->dxyz.x = pj_param (P->ctx, P->params, "ddx").f;

    if (pj_param(P->ctx, P->params, "tdy").i)
        Q->dxyz.y = pj_param (P->ctx, P->params, "ddy").f;

    if (pj_param(P->ctx, P->params, "tdz").i)
        Q->dxyz.z = pj_param (P->ctx, P->params, "ddz").f;

    /* Rotations rates */
    if (pj_param (P->ctx, P->params, "tdrx").i)
        Q->dopk.o = pj_param (P->ctx, P->params, "ddrx").f * ARCSEC_TO_RAD;

    if (pj_param (P->ctx, P->params, "tdry").i)
        Q->dopk.p = pj_param (P->ctx, P->params, "ddry").f * ARCSEC_TO_RAD;

    if (pj_param (P->ctx, P->params, "tdrz").i)
        Q->dopk.k = pj_param (P->ctx, P->params, "ddrz").f * ARCSEC_TO_RAD;

    if (pj_param (P->ctx, P->params, "tdtheta").i)
        Q->dtheta = pj_param (P->ctx, P->params, "ddtheta").f * ARCSEC_TO_RAD;

    /* Scale rate */
    if (pj_param (P->ctx, P->params, "tds").i)
        Q->dscale = pj_param (P->ctx, P->params, "dds").f;


    /* Epoch */
    if (pj_param(P->ctx, P->params, "tt_epoch").i)
        Q->t_epoch = pj_param (P->ctx, P->params, "dt_epoch").f;

    Q->xyz    =  Q->xyz_0;
    Q->opk    =  Q->opk_0;
    Q->scale  =  Q->scale_0;
    Q->theta  =  Q->theta_0;

    if ((Q->opk.o==0)  && (Q->opk.p==0)  && (Q->opk.k==0) && (Q->scale==0) &&
        (Q->dopk.o==0) && (Q->dopk.p==0) && (Q->dopk.k==0)) {
        Q->no_rotation = 1;
    }

    if( !read_convention(P) ) {
        return nullptr;
    }

    /* Let's help with debugging */
    if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_DEBUG) {
        proj_log_debug(P, "Helmert parameters:");
        proj_log_debug(P, "x=  %8.5f  y=  %8.5f  z=  %8.5f", Q->xyz.x, Q->xyz.y, Q->xyz.z);
        proj_log_debug(P, "rx= %8.5f  ry= %8.5f  rz= %8.5f",
                Q->opk.o / ARCSEC_TO_RAD, Q->opk.p / ARCSEC_TO_RAD, Q->opk.k / ARCSEC_TO_RAD);
        proj_log_debug(P, "s=  %8.5f  exact=%d%s", Q->scale, Q->exact,
                       Q->no_rotation ? "" :
                       Q->is_position_vector ? "  convention=position_vector" :
                       "  convention=coordinate_frame");
        proj_log_debug(P, "dx= %8.5f  dy= %8.5f  dz= %8.5f",   Q->dxyz.x, Q->dxyz.y, Q->dxyz.z);
        proj_log_debug(P, "drx=%8.5f  dry=%8.5f  drz=%8.5f",   Q->dopk.o, Q->dopk.p, Q->dopk.k);
        proj_log_debug(P, "ds= %8.5f  t_epoch=%8.5f", Q->dscale, Q->t_epoch);
    }

    if (Q->no_rotation) {
        return P;
    }

    update_parameters(P);
    build_rot_matrix(P);

    return P;
}


/***********************************************************************/
PJ *TRANSFORMATION(molobadekas, 0) {
/***********************************************************************/

    struct pj_opaque_helmert *Q;

    if( !init_helmert_six_parameters(P) ) {
        return nullptr;
    }

    P->fwd3d  = helmert_forward_3d;
    P->inv3d  = helmert_reverse_3d;

    Q = (struct pj_opaque_helmert *)P->opaque;

    /* Scale */
    if (pj_param (P->ctx, P->params, "ts").i) {
        Q->scale_0 = pj_param (P->ctx, P->params, "ds").f;
    }

    Q->opk    =  Q->opk_0;
    Q->scale  =  Q->scale_0;

    if( !read_convention(P) ) {
        return nullptr;
    }

    /* Reference point */
    if (pj_param (P->ctx, P->params, "tpx").i)
        Q->refp.x = pj_param (P->ctx, P->params, "dpx").f;

    if (pj_param (P->ctx, P->params, "tpy").i)
        Q->refp.y = pj_param (P->ctx, P->params, "dpy").f;

    if (pj_param (P->ctx, P->params, "tpz").i)
        Q->refp.z = pj_param (P->ctx, P->params, "dpz").f;


    /* Let's help with debugging */
    if (proj_log_level(P->ctx, PJ_LOG_TELL) >= PJ_LOG_DEBUG) {
        proj_log_debug(P, "Molodensky-Badekas parameters:");
        proj_log_debug(P, "x=  %8.5f  y=  %8.5f  z=  %8.5f", Q->xyz_0.x, Q->xyz_0.y, Q->xyz_0.z);
        proj_log_debug(P, "rx= %8.5f  ry= %8.5f  rz= %8.5f",
                Q->opk.o / ARCSEC_TO_RAD, Q->opk.p / ARCSEC_TO_RAD, Q->opk.k / ARCSEC_TO_RAD);
        proj_log_debug(P, "s=  %8.5f  exact=%d%s", Q->scale, Q->exact,
                       Q->is_position_vector ? "  convention=position_vector" :
                       "  convention=coordinate_frame");
        proj_log_debug(P, "px= %8.5f  py= %8.5f  pz= %8.5f",   Q->refp.x, Q->refp.y, Q->refp.z);
    }

    /* as an optimization, we incorporate the refp in the translation terms */
    Q->xyz_0.x +=  Q->refp.x;
    Q->xyz_0.y +=  Q->refp.y;
    Q->xyz_0.z +=  Q->refp.z;

    Q->xyz    =  Q->xyz_0;

    build_rot_matrix(P);

    return P;
}

template<class UnaryFunction>
void recursive_iterate(const json& j, UnaryFunction f)
{
	for (auto it = j.begin(); it != j.end(); ++it)
	{
		// auto typenam = it->type_name();
		auto v = it.value();

		if (it->is_array() || it->is_object())
		{
			recursive_iterate(*it, f);
		}	
		else if (it->is_null())
		{
			f(it);
		}
		else if (it->is_number_float())
		{
			auto value = it.value();
		}		 
		else
		{
			f(it);
		}
	}
}

// TODO: Flytte all GeoJson til eigen klasse.
/***********************************************************************
*
/***********************************************************************/
static void testReadGeojson(/*char* fileName*/)
{
	char* fileName = "C:/Prosjekter/SkTrans/EurefNgo/Punksky_tilfeldig/Ngo_areas.geojson";
 
	FILE *f = fopen(fileName, "rb");
	fseek(f, 0, SEEK_END);
	long fsize = ftell(f);
	fseek(f, 0, SEEK_SET);

	char *string = (char *)malloc(fsize + 1);
	fread(string, fsize, 1, f);
	fclose(f);

	string[fsize] = 0;
	 
	// parse and serialize JSON
	json j_complete = json::parse(string);
	// std::cout << std::setw(4) << j_complete << "\n\n";
		 
	recursive_iterate(j_complete, [](json::const_iterator it)
	{});

	// Test
	auto feat = j_complete.at("features");
	//auto c = j_complete.find("coordinates");
 
	for (auto it = feat.begin(); it != feat.end(); ++it)
	{
		//auto hh =	it1.key["geometry"];

		auto geo = (*it)["geometry"];
		auto coords = geo.at("coordinates");

		//recursive_iterateToFloat(*it, coords.value);

		//if (it1->is_number_float)
		/*
		for (auto it2 = geo.begin(); it2 != geo.end(); ++it2)
		{
			auto coord = (*it2)["coordinates"];
		}*/
	}

	for (auto& x : feat.items())
	{
		//std::cout << "key: " << x.key() << ", value: " << x.value() << '\n';
	}	 
	/*
	 json::iterator it = feat.begin();

	 std::for_each(feat.begin(), feat.end(), [](std::string &string)
	 { 
		 std::cout << string << "\n";
	 }
	 );
*/
	//feat.at()
	auto feat1 = feat.at(0);
	auto geo = feat1.at("geometry");
	auto coords1 = geo.at("coordinates");
	auto coords2 = coords1.at(0);
	auto coords3 = coords2.at(0);
	auto coords4 = coords3.at(0);
	auto coords5 = coords4.at(0);

	//auto multiPolygon = j_complete.array("MultiPolygon");
	//auto items = j_complete.items();
	auto json_string = j_complete.dump();
}

/******************************************************************************************
* http://www.mygeodesy.id.au/documents/Coord%20Transforms%20in%20Cadastral%20Surveying.pdf
*
* https://www.degruyter.com/downloadpdf/j/rgg.2014.97.issue-1/rgg-2014-0009/rgg-2014-0009.pdf
/******************************************************************************************/
static void calculateHelmertParameters(std::vector<CommonPointPair> *commonPointList, PJ_LP lp)
{
	auto n = commonPointList->size();
	 
	double k = 0.00039;
	double c = 0.06900;

	double coslat = cos(lp.phi * M_PI / 180.0);

	double x = lp.phi;
	double y = lp.lam * coslat;		 

	MatrixXd cnn(n, n);
	MatrixXd cmn(n, 1);

	// Vector From System
	MatrixXd xF(n, 1);
	MatrixXd yF(n, 1);

	// Vector To System
	MatrixXd xT(n, 1);
	MatrixXd yT(n, 1);	
	
	for (int i = 0; i < n; i++)
	{
		CommonPointPair point1 = commonPointList->at(i);
		xF(i, 0) = point1.fromPoint.phi;
		yF(i, 0) = point1.fromPoint.lam * coslat;

		xT(i, 0) = point1.toPoint.phi;
		yT(i, 0) = point1.toPoint.lam * coslat;		
	}
    // cout << "xF:" << endl << xF << endl;
 	//cout << "yF:" << endl << yF << endl;
	//cout << "xT:" << endl << xT << endl;
	//cout << "yT:" << endl << yT << endl;

	for (int i = 0; i < n; i++)
	{
		double dist = hypot(xF(i, 0) - x, yF(i, 0) - y);
		double a = (M_PI / 2.0) * (dist / c);
		cmn(i, 0) = k * exp(-a) * cos(a);

		for (int j = 0; j < n; j++)
		{
			dist = hypot(xF(i, 0) - xF(j, 0), yF(i, 0) - yF(j, 0));
			a = (M_PI / 2.0) * (dist / c);
			cnn(i, j) = k * exp(-a) * cos(a);
		}
	}
	//	cout << "cmn:" << endl << cmn << endl;
	//	cout << "cnn:" << endl << cnn << endl;

	MatrixXd p = cnn.inverse();
	//cout << "p:" << endl << p << endl;
	 
	// N, diagonal sum of p:
	MatrixXd sep = p * MatrixXd::Ones(n, 1);
	MatrixXd N = MatrixXd::Ones(1, n) * sep;
	//cout << "N:" << endl << N << endl;

	// Mean values
	MatrixXd xFT = MatrixXd::Ones(1, n) * p * xF;
	MatrixXd yFT = MatrixXd::Ones(1, n) * p * yF;
	MatrixXd xTT = MatrixXd::Ones(1, n) * p * xT;
	MatrixXd yTT = MatrixXd::Ones(1, n) * p * yT;
	 
	// Mass center: 
	MatrixXd xF0 = xFT * N.inverse();
	MatrixXd yF0 = yFT * N.inverse();
	MatrixXd xT0 = xTT * N.inverse();
	MatrixXd yT0 = yTT * N.inverse();
	
	// Coordinates with Mass center origin
	MatrixXd dxF = xF - MatrixXd::Ones(n, 1) * xF0;
	MatrixXd dyF = yF - MatrixXd::Ones(n, 1) * yF0;
	MatrixXd dxT = xT - MatrixXd::Ones(n, 1) * xT0;
	MatrixXd dyT = yT - MatrixXd::Ones(n, 1) * yT0;

	double n11 = 0.0;
	double t1 = 0.0;
	double t2 = 0.0;

	for (int i = 0; i < n; i++)
	{
		n11 += (pow(dxF(i), 2) * sep(i)) + (pow(dyF(i), 2) * sep(i));
		t1 += (dxF(i) * dxT(i) + dyF(i) * dyT(i)) * sep(i);
		t2 += (dyF(i) * dxT(i) - dxF(i) * dyT(i)) * sep(i);
	}
	double a = t1 / n11;
	double b = t2 / n11;
 
	double tx = xT0(0) - a * xF0(0) - b * yF0(0);
	double ty = yT0(0) + b * xF0(0) - a * yF0(0);

	// Sigma noise
	MatrixXd snx(n, 1);
	MatrixXd sny(n, 1);

	for (int i = 0; i < n; i++)
	{
		snx(i) = dxT(i) - a * dxF(i) - b * dyF(i);
		sny(i) = dyT(i) + b * dxF(i) - a * dyF(i);
	}
	//cout << "snx:" << endl << snx << endl;
	//cout << "sny:" << endl << sny << endl;

	double xTrans = xT0(0) - a * (xF0(0) - x) - b * (yF0(0) - y);
	double yTrans = yT0(0) + b * (xF0(0) - x) - a * (yF0(0) - y);
	
	MatrixXd smx = cmn.transpose() * p * snx;
	MatrixXd smy = cmn.transpose() * p * sny;
	//cout << "smx:" << endl << smx << endl;
	//cout << "smy:" << endl << smy << endl;

	double xEst = xTrans + smx(0);
	double yEst = yTrans + smy(0);
}
 
bool DistanceLess(const CommonPointPair& lhs, const CommonPointPair& rhs)
{
	return lhs.dist < rhs.dist;
}

/***********************************************************************
* https://stackoverflow.com/questions/4509798/finding-nearest-point-in-an-efficient-way
/***********************************************************************/
std::vector<CommonPointPair> findClosestPoints(std::vector<CommonPointPair> *commonPointList, PJ_LP lp, int n, int areaId)
{
	std::vector<CommonPointPair> distances;
	std::vector<CommonPointPair> closestDistances;
	double coslat = cos(lp.phi * M_PI / 180.0);

	for each (CommonPointPair pair in *(commonPointList))
	{
		double deltaPhi = pair.fromPoint.phi - lp.phi;
		double deltaLam = (pair.fromPoint.lam - lp.lam) * coslat;
		
		pair.dist = sqrt ((deltaPhi * deltaPhi) + (deltaLam * deltaLam));
		distances.push_back(pair);
	}
	
	std::sort(distances.begin(), distances.end(), DistanceLess);

	for (int i = 0; i < n; i++)
		closestDistances.push_back(distances[i]);

	return closestDistances;
}

/***********************************************************************
*
/***********************************************************************/
bool PointIsInArea(PJ_LP pointPJ_LP, char* fileName)
{	 
	std::ifstream file(fileName, std::ios::in);

	Point points[] = { {} };
	Point point = { pointPJ_LP.phi, pointPJ_LP.lam };
	vector<Point> pointVector;

	if (file.is_open())
	{
		std::vector<std::vector<std::string> > dataList;
		std::string line = "";
		double x, y;

		while (!file.eof())
		{
			getline(file, line, ',');
			x = atof(line.c_str());

			getline(file, line, '\n');
			y = atof(line.c_str());

			Point areaPoint;
			areaPoint.x = x;
			areaPoint.y = y;

			pointVector.push_back(areaPoint);
		}
		file.close();
	};

	Point *vectorPointer = pointVector.data();
	int n = size(pointVector);

	return isInside(vectorPointer, n, point);
}

/***********************************************************************
*
/***********************************************************************/
int AreaIdPoint(PJ_LP pointPJ_LP) // TODO: Endre namn og argument 
{	
	// TODO: Flytte områdefilene
	char* fileName2 = "C:/Prosjekter/SkTrans/EurefNgo/Punksky_tilfeldig/Area2.csv";
	char* fileName3 = "C:/Prosjekter/SkTrans/EurefNgo/Punksky_tilfeldig/Area3.csv";
	char* fileName4 = "C:/Prosjekter/SkTrans/EurefNgo/Punksky_tilfeldig/Area4.csv";
	
	// TODO: Områdefil som geojson. Json ligg under include/proj/internal/nlohmann
	//testReadGeojson();

	if (PointIsInArea(pointPJ_LP, fileName2))
		return 2;
	else if (PointIsInArea(pointPJ_LP, fileName3))
		return 3;
	else if (PointIsInArea(pointPJ_LP, fileName4))
		return 4;

	return 1;
}

/***********************************************************************
*
/***********************************************************************/
PJ_LP proj_commonPointInit(PJ_LP lp)
{
	std::vector<CommonPointPair> commonPointList;

	// TODO: Flytte lan1_fellesp_20081014.cpt til ei anna mappe
	char* fileName = "C:/Users/Administrator/source/repos/Skproj/Octave/lan1_fellesp_20081014.cpt";
	std::ifstream file(fileName, std::ios::in | std::ios::binary | std::ios::ate);

	int areaId = AreaIdPoint(lp);

	if (file.is_open())
	{
		int bufferSize8 = 8; 
		int bufferSize4 = 4;
		std::string name = "";
		char *charBuffer8 = new char[bufferSize8];
		char *charBuffer4 = new char[bufferSize4];	

		file.seekg(0, std::ios::beg);
		
		while (!file.eof())
		{
			CommonPointPair	commonPoint;

			file.read(charBuffer8, bufferSize8);
			name = "";

			for (int i = 0; i < bufferSize8; i++) 			
				name += charBuffer8[i];

			commonPoint.name = name;
			
			file.read(charBuffer8, bufferSize8);
			commonPoint.fromPoint.phi = *(reinterpret_cast<double*>(charBuffer8));

			file.read(charBuffer8, bufferSize8);
			commonPoint.fromPoint.lam = *(reinterpret_cast<double*>(charBuffer8));

			file.read(charBuffer8, bufferSize8);
			commonPoint.toPoint.phi = *(reinterpret_cast<double*>(charBuffer8));

			file.read(charBuffer8, bufferSize8);
			commonPoint.toPoint.lam = *(reinterpret_cast<double*>(charBuffer8));

			file.read(charBuffer4, bufferSize4);
			commonPoint.area = *(reinterpret_cast<__int32*>(charBuffer4));

			if (areaId != commonPoint.area)
				continue;

			commonPointList.push_back(commonPoint);
		}
		file.close();
	};

    int numberOfSelectedPoints = 20; 

	// TODO: Flytte kalla	
	auto closestPoints = findClosestPoints(&commonPointList, lp, numberOfSelectedPoints, areaId);

	// TODO: Leggje inn Helmert her.
	calculateHelmertParameters(&closestPoints, lp);

	return lp;
}
