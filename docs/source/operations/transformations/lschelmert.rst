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

    A link to geojson multipolygons. The operation tests if the input coordinates are
	within some of the multipolygons. Multipolygons have a foreignkey areaid which
	is a field in the point pair object from the cpt-file. Point pairs are selected 
	based on selected multipolygon.

.. option:: +points=<value>

    The number of maximum selected point candidates used in Least Square Collocation 
	and 2D Helmert. Default is 20.  Units of latitude and longitude is in radians, and 
	height in meters.

.. option:: +maximum_dist=<value>

    The maximum distance between input point and selected point candidate. Unit of the
	distance is radians. Default is 0.1 radians.

.. option:: +c_coll=<value>

    ....
	Default is 0.00039.

.. option:: +k_coll=<value>
    
	....
	Default is 0.001204.
 

Mathematical description
################################################################################


Least Squared Colloction
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


2D Helmert transformation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

