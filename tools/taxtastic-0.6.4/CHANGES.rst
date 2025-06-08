==========================
 change log for taxtastic
==========================

0.6.4
=====
  * rank cohort has been added as an official rank
  * ``taxit taxtable --taxtable`` has been added back to work from a pre-built taxtable.  This switch used
    to be called ``--from-table``

0.6.3
=====
  * root is now a valid node

0.6.2
=====
  * bug fixes

0.6.1
=====
  * bug fixes

0.6.0
=========
 * ``taxit update_taxids --taxid-column`` allows updating of any tax_id column [GH-84]
 * ``taxit update_taxids --ignore-unknowns`` allows unknown tax_ids to remain in final output [GH-84]
 * ``taxit new_database`` adds ncbi as the only entry in the source table [GH-91]
 * ``taxit new_database`` adds a ranks table with all ranks appearing in the nodes table [GH-86]
 * Any flavor of database can used with Taxtastic.  But only sqlite and Postgres have been tested.
 * Numerous new features and performance improvements.

0.5.7
=====
 * ``taxit update_taxids`` is significantly faster but can still use some optimizations [GH-78]

0.5.6
=====
 * remove support for reading excel spreadsheets (GH-71)
 * requirements.txt identifies all direct dependencies
 * add ``taxit merge_taxtables``
 * New ``taxit new_database --taxdump`` and ``taxit new_database--taxdump-url`` arguments
   for flexibility on taxdump.zip location(s)
 * New function taxtable.remove_subtree() (GH-80)

0.5.5
=====
 * new ``taxit count_taxids`` counts every tax_id occurance in a ``taxit taxtable`` lineage [GH-75]
 * new ``taxit taxid_classified`` decides if a tax_id is primary and valid (True/False)
 * ``taxit update_taxids`` will halt on unknown tax_ids unless ``--unknowns FILE`` is specified
 * ``taxit update_taxids`` only requires a csv file with 'tax_id' column
 * ``taxit update_taxids`` takes an optional ``--name-column`` to assist in assigning tax_ids
 * ``taxit update_taxids`` will read stdin if csv file is not provided as argument

0.5.4
=====

 * Add ``taxit taxtable --full`` for outputing all ranks in header.
 * Update subcommand help text
 * Generate Sphinx docs using help text emitted by subcommands (GH-70)

0.5.3
=====

 * Suppress warning when updating refpkg ``tree_stats`` file via ``taxit update``.

0.5.2
=====

 * Fix GH-63: "empirical_frequencies" now set to false when parsing FastTree AA statistics files
 * Close GH-64: "empirical_frequencies" is now available as a flag for PhyML statistics files
 * Fixed bug that prevented temporary files from being deleted

0.5.1
=====

 * Fix GH-62: "empirical_frequencies" was not set when parsing PhyML AA statistics files.

0.5.0
=====

 * Add ``.drop()`` ``.collapse()`` methods to ``taxtastic.taxtable.TaxNode``
 * Change ``is_classified`` column in taxonomy database: now does not mark
   below species as unclassified if the species-level classification is valid. [GH-59]
 * Add ``taxit composition`` - shows the taxonomic composition of a reference package at a given rank
 * Fix broken ``taxit lonelynodes``
 * Add ``taxit merge`` - Identifies tax_ids which have been merged, suggests new tax_ids.
 * Add ``taxit add_to_taxtable`` - adds nodes to a taxonomy [GH-60]
 * Fix support for newer versions of PhyML [GH-61]
 * Updates for compatibility with RAxML 7.7.2


0.4
===

 * 'names' table in the taxonomy database has a new column
   'is_classified' indicating whether 'tax_name' should be considered
   "classified".
 * Bugfix in ``taxit findcompany``
 * Support stdin as a source for ``taxit findcompany``
 * Taxonomy objects use NCBI ranks by default
 * Reference packages are created optionally (fixes creation of empty reference
   packages for commands like ``taxit info nonexisting.refpkg``)
 * Support zipped reference packages
 * Add a taxtable API: ``taxtastic.taxtable``
 * Remove some tests requiring a full taxonomy database
 * Rerooting reference packages on creation [GH-57]
 * More intelligent file name generation on clash
 * Deprecate the default ``create=True`` in ``taxtastic.refpkg.Refpkg``
 * Some PEP8 fixes


0.3.2
=====

 * version number contains abbreviated git sha identifying the commit.
 * Initial release to PyPI
 * Add findcompany subcommand
 * Add refpkg_intersection subcommand
 * Remove some obsolete components
 * Check required fields in seqinfo file [GH-46]
 * Add option to build taxtable from seqinfo file [GH-55]
 * Add subcommand to update taxids [GH-56]
 * Support FastTree AA and DNA log files
 * Fix rank order bug (infraorder was below parvorder)
 * Documentation updates
