.. _lschelmert:

================================================================================
Least Squares Collaction with Helmert 2D Transformation
================================================================================

.. versionadded:: 7.1.0

.. note::
	The operation Helmert 2D Transformation Least with Squares Collocation is
	processes once for each point. 

+---------------------+----------------------------------------------------------+
| **Alias**           | lschelmert                                               |
+---------------------+----------------------------------------------------------+
| **Domain**          | 2D                                                       |
+---------------------+----------------------------------------------------------+
| **Input type**      | Geodetic coordinates                                     |
+---------------------+----------------------------------------------------------+
| **output type**     | Geodetic coordinates                                     |
+---------------------+----------------------------------------------------------+


Mathematical description
################################################################################

When the Norwegian Mapping Authority introduced the new geodetic datum EUREF89,
it was necessary to find a proper transformation technic beetween  the deprecated
geodetic datum NGO1948 to the EUREF89. NGO1948 was significant deformated, hereby 
it was computed and realized in different areas of the country. A consequence of
this is huge planar gaps in some counties and municipalities. The gaps are upon
2-3 meters on borders. Description and evaluation of the method are further
documented in the articles see :cite:`OMathisen2002` and :cite:`OMathisen2003`.

The processing is done in three steps:

	1. Selection of point pairs
	2. Deterministic step: 2D Helmert transformation
	3. Statistic step: Smoothing Least Squared Collocation

Selection of point pairs
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Both computations in step 2 and 3 are based on a certain number of point pairs
according to the parameter ``+points``. The method selects the ``+points`` number
of points that are closest to the input point. The selected points have to be
within a radius of ``+maximum_dist`` km.


2D Helmert transformation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In the first step 2D Helmert transformation parameters are computed based on a
certain number of selected common points (point pairs). 2D Helmert transformation
consist of four parameters and those are solved by traditional Least Squares Method.
The inverted covariance matrix from LSC is adopted as weight matrix. That means
closer points gets higher weight.
 
A standard 2D Helmert is described as:

.. math::
    :label: helmert2d

    \left[\begin{array}{cc}
    x \\
    y
    \end{array}\right]+\left[\begin{array}{cc}
    v_x \\
    v_y
    \end{array}\right]=\left[\begin{array}{cc}
    a & b \\
    -b & a
    \end{array}\right]\left[\begin{array}{cc}
    u \\
    v
    \end{array}\right]+\left[\begin{array}{cc}
    t_x \\
    t_y
    \end{array}\right]\\

Where :math:`u` og :math:`v` is 2D coordinates in source coodinate system and x og y in
target coordinate system. Helmert transformation parametres are denoted as
:math:`t_x`, :math:`t_y`, :math:`a` and :math:`b`.

The selected covariance function for this operation a modified first Gauss Markov.


Covariance matrix of the given common points:

.. math::
    :label: cov_nn

    C_{nn}=ke^{-\frac{\pi{}}{2}\frac{d}{c}}\cos{\frac{\pi{}}{2}\frac{d}{c}}\\*


where :math:`n` is the number of common points, :math:`d` distance in km, :math:`c` the ccoll parameter and :math:`k` is the kcoll parameter.


Covariance matrix of the input point:

.. math::
    :label: cov_mn

    C_{mn}=ke^{-\frac{\pi{}}{2}\frac{d}{c}}\cos{\frac{\pi{}}{2}\frac{d}{c}}\\*


where :math:`n` is the number of common points, :math:`m` is the number of transformed and predicted points. :math:`m` is basically. :math:`d` is distance in km, :math:`c` the ccoll parameter and :math:`k` the kcoll parameter.


Further mass center points are computed for both coordinate systems with weight from the inverted covariance function. The weights are noted :math:`W`.


Weight matrix, inverse of Cnn:

.. math::
    :label: weight_mat

    W={C_{nn}}^{-1}

\

Ws is the sum of the entired weight matrix:


.. math::

    w_s=\sum_{i=1}^n\sum_{j=1}^nw_{ji}

\

Sum weight for each point:

.. math::

    w=W\ \vec{1}

\
 
Mass center computed based on weighed centroid:

.. math::

    \begin{array}{cc}u_0=\frac{w^Tu}{w_s}\end{array}
    \begin{array}{cc}v_0=\frac{w^Tv}{w_s}\end{array}
    \begin{array}{cc}x_0=\frac{w^Tx}{w_s}\end{array}
    \begin{array}{cc}y_0=\frac{w^Ty}{w_s}\end{array}

\

Target and source points moved to mass center as centroids:

.. math::

    \begin{array}{cc}\bar{u}=u-\vec{1}u_0\end{array}
    \begin{array}{cc}\bar{v}=v-\vec{1}v_0\end{array}
    \begin{array}{cc}\bar{x}=x-\vec{1}x_0\end{array}
    \begin{array}{cc}\bar{y}=y-\vec{1}y_0\end{array}

\

The modified observation equation is now transformed with centroids as input and output:

.. math::
    :label: helmert2d_mod

    \left[\begin{array}{cc}
    \bar{x} \\
    \bar{y}
    \end{array}\right]+\ \left[\begin{array}{cc}
    v_x \\
    v_y
    \end{array}\right]=\left[\begin{array}{cc}
    a & b \\
    -b & a
    \end{array}\right]\left[\begin{array}{cc}
    \bar{u} \\
    \bar{v}
    \end{array}\right]+\left[\begin{array}{cc}
    T_x \\
    T_y
    \end{array}\right]

\

Least Squares Estimation of Helmert 2D parameter based on simplified inversed normal equation:

.. math::
    :label: normal_eq

    \left[\begin{array}{cc}
    \sum_{i=1}^nw_i({{\bar{u}}_i}^2+{{\bar{v}}_i}^2) & 0 \\
    0 & \sum_{i=1}^nw_i({{\bar{u}}_i}^2+{{\bar{v}}_i}^2)
    \end{array}\right]\left[\begin{array}{cc}
    a \\
    b
    \end{array}\right]=\left[\begin{array}{cc}
    \sum_{i=1}^nw_i({\bar{u}}_i{\bar{x}}_i+{\bar{v}}_i{\bar{y}}_i) \\
    \sum_{i=1}^nw_i({\bar{v}}_i{\bar{x}}_i-{\bar{u}}_i{\bar{y}}_i)
    \end{array}\right]

\

Solved Helmert scale/rotation parameters :math:`a` and :math:`b`:

.. math::
    :label: normal_ab

    \begin{array}{cc}a=\frac{\sum_{i=1}^nw_i({\bar{u}}_i{\bar{x}}_i+{\bar{v}}_i{\bar{y}}_i)}{\sum_{i=1}^nw_i({{\bar{u}}_i}^2+{{\bar{v}}_i}^2)}\end{array}
    \begin{array}{cc}b=\frac{\sum_{i=1}^nw_i({\bar{v}}_i{\bar{x}}_i-{\bar{u}}_i{\bar{y}}_i)}{\sum_{i=1}^nw_i({{\bar{u}}_i}^2+{{\bar{v}}_i}^2)}\end{array}

\

Solving Helmert translation parameters :math:`t_x`, :math:`t_y`:

.. math::
    :label: normal_t

    \begin{array}{cc}t_x=x_0-u_0a-v_0b\end{array}
    \begin{array}{cc}t_y=y_0+u_0b-v_0a\end{array}

\

Residuals from least squares 2D Helmert:

.. math::
    :label: residual_xy

    \begin{array}{cc}v_x=\bar{x}-a\bar{u}-b\bar{v}\end{array}
    \begin{array}{cc}v_y=\bar{y}+b\bar{u}-a\bar{v}\end{array}

\

Input coordinate transformed to the target coordinate system:

.. math::
    :label: pred_xy

    {\varphi{}}_H=x_0-a\left(u_0-{\varphi{}}_{in}\right)-b(v_0-{\lambda{}}_{in}\cos{{\varphi{}}_{in}})

    {\lambda{}}_H=\frac{y_0+b\left(u_0-{\varphi{}}_{in}\right)-a(v_0-{\lambda{}}_{in}\cos{{\varphi{}}_{in}})}{\cos{{\varphi{}}_{in}}}
 

Least Squared Collocation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The signal of the given common points are set as the same as the computed residuals from
the least squares 2D Helmert.

.. math::
    :label: signal_n

    \begin{array}{cc}s_{nx}=v_x\end{array}
    \begin{array}{cc}s_{ny}=v_y\end{array}

\

Then the signal of the transformed points is given by:

.. math::
    :label: signal_m

    \begin{array}{cc}s_{mx}=C_{mn}W\ s_{nx}\end{array}
    \begin{array}{cc}s_{my}=C_{mn}W\ s_{ny}\end{array}

\
 
The signal from Least Squares Collocation is added to the tranformed point. The location is called predicted point.

\

Predicted output latitude:

.. math::
    :label: phi_out

    {\varphi{}}_{out}={\varphi{}}_H+s_{mx}

\

Predicted output longitude:

.. math::
    :label: lambda_out

    {\lambda{}}_{out}={\lambda{}}_H+\frac{s_{my}}{\cos{{\varphi{}}_{in}}}


Examples
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The most simple proj string refers to a binary file with list of common points.
A common point is defined by two sets of coordinates, one for the source 
coordinate system and one for the target coodinate system:

::

    +proj=lschelmert
    +pp_trans=EUREF89_NGO48_20081014.cpt


By adding the parameter `+polygons`, the selection of points might be separated in different areas:

::

    +proj=lschelmert
    +pp_trans=EUREF89_NGO48_20081014.cpt
    +polygons=Flater.geojson
    +ellps=GRS80

Proj string with entired set of optional parameters:

::

    +proj=lschelmert
    +pp_trans=EUREF89_NGO48_20081014.cpt
    +polygons=Flater.geojson
    +points=15
    +maximum_dist=80.0
    +ccoll=10.0
    +kcoll=0.0005

Parameters
###############################################################################

Required
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. option:: +pp_trans=<list>

    A link to file with list of point pairs. A point pair is a object with coordinates referred in two geodetic datums.
    The file itselfs is in binary format.

    If a file is prefixed by an @ the file is considered optional and PROJ will the not complain if the file is not available.

Optional
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. option:: +polygons=<list>

    A link to geojson multipolygons. The operation tests if the input coordinates are within some of the multipolygons.
    Multipolygons have a foreignkey areaid which is a field in the point pair object from the cpt-file.
    Point pairs are selected based on selected multipolygon.

    If a file is prefixed by an @ the file is considered optional and PROJ will the not complain if the file is not available.

.. option:: +points=<value>

    The number of maximum selected point candidates used in Least Square Collocation and 2D Helmert.
    Units of latitude and longitude is in radians, and height in meters.

    *Default is 20.*

.. option:: +maximum_dist=<value>

    The maximum distance between input point and selected point candidate. Unit of the distance is km. 

    *Default is 100.0 km.*

.. option:: +ccoll=<value>
    
    The ccoll value is the distance where the empirical covariance touches zero. The unit ccoll is in km. 

    *Default is 7.7.*

.. option:: +kcoll=<value>

    The kcoll coefficient is simular to C0 in a standard Gauss Markov first order covariance function.

    *Default is 0.00039.*

.. include:: ../options/ellps.rst


Further reading
###############

#. `Ligas, Banasik "LSC Alternative to Helmert's Tranformation with Hausbrandt's Post-Tranformation Correction" <https://www.degruyter.com/downloadpdf/j/rgg.2014.97.issue-1/rgg-2014-0009/rgg-2014-0009.pdf>`_

#. `R. E. Deakin "Coordinate Transformation for Cadastral Surveying" <http://www.mygeodesy.id.au/documents/Coord%20Transforms%20in%20Cadastral%20Surveying.pdf>`_

#. `Stackoverflow.com <https://stackoverflow.com/questions/4509798/finding-nearest-point-in-an-efficient-way>`_

