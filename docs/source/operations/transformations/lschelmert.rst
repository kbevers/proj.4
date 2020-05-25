.. _lschelmert:

================================================================================
Least Squares Collaction with Helmert 2D Transformation
================================================================================

.. versionadded:: X.X.X

.. note::


+---------------------+----------------------------------------------------------+
| **Alias**           | lschelmert                                               |
+---------------------+----------------------------------------------------------+
| **Domain**          | 2D                                                       |
+---------------------+----------------------------------------------------------+
| **Input type**      | Geodetic coordinates                                     |
+---------------------+----------------------------------------------------------+
| **output type**     | Geodetic coordinates                                     |
+---------------------+----------------------------------------------------------+

Examples
###############################################################################

    +proj=lschelmert +pp_trans=EUREF89_NGO48_20081014.cpt +polygons=Flater.geojson +ellps=GRS80


Parameters
################################################################################

Required
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. option:: +pp_trans=<File>

    A link to file with list of point pairs. A point pair is a object with coordinates
	referred in two geodetic datums. The file itselfs is in binary format.

Optional
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. option:: +polygons=<File>

    A link to geojson multipolygons. The operation tests if the input coordinates
	are within some of the multipolygons. Multipolygons have a foreignkey areaid
	which is a field in the point pair object from the cpt-file. Point pairs are
	selected based on selected multipolygon.

.. option:: +points=<value>

    The number of maximum selected point candidates used in Least Square Collocation
	and 2D Helmert. Default is 20.  Units of latitude and longitude is in radians, and
	height in meters.

.. option:: +maximum_dist=<value>

    The maximum distance between input point and selected point candidate. Unit of the
	distance is radians. Default is 0.1 radians.

.. option:: +c_coll=<value>
    
	The c_coll value is the distance where the empirical covariance touches zero. The
	unit c_coll is in km.

    Default is 7.7.
	

.. option:: +k_coll=<value>

    The k_coll coefficient is simular to C0 in a standard Gauss Markov first order covariance
	function.
	
	Default is 0.00039.


Mathematical description
################################################################################

When the Norwegian Mapping Authority introduced the new geodetic datum EUREF89,
it was necessary to find a proper transformation technic beetween  the deprecated
geodetic datum NGO1948 to the EUREF89. NGO1948 was significant deformated, hereby 
it was computed and realized in different areas of the country. A consequence of
this is huge planar gaps in some counties and municipalities. The gaps are upon
2-3 meters on borders. Description and evaluation of the method are further
documented in the articles see :cite:`OMathisen2002` and :cite:`OMathisen2003`.

2D Helmert transformation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In the first step 2D Helmert transformation parameters are computed based on a
certain number of selected common points (point pairs).




Least Squared Colloction
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The selected covariance function for this operation a modified first Gauss Markov.

.. math::
    :label: covfunc
	
	C{d} = X_{in} + \Delta X
	 


