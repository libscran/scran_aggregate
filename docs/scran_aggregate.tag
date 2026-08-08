<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.12.0">
  <compound kind="file">
    <name>aggregate_across_cells.hpp</name>
    <path>scran_aggregate/</path>
    <filename>aggregate__across__cells_8hpp.html</filename>
    <class kind="struct">scran_aggregate::AggregateAcrossCellsOptions</class>
    <class kind="struct">scran_aggregate::AggregateAcrossCellsBuffers</class>
    <class kind="struct">scran_aggregate::AggregateAcrossCellsResults</class>
    <namespace>scran_aggregate</namespace>
  </compound>
  <compound kind="file">
    <name>aggregate_across_genes.hpp</name>
    <path>scran_aggregate/</path>
    <filename>aggregate__across__genes_8hpp.html</filename>
    <class kind="struct">scran_aggregate::AggregateAcrossGenesOptions</class>
    <class kind="struct">scran_aggregate::AggregateAcrossGenesBuffers</class>
    <class kind="struct">scran_aggregate::AggregateAcrossGenesSet</class>
    <class kind="struct">scran_aggregate::AggregateAcrossGenesResults</class>
    <namespace>scran_aggregate</namespace>
  </compound>
  <compound kind="file">
    <name>scran_aggregate.hpp</name>
    <path>scran_aggregate/</path>
    <filename>scran__aggregate_8hpp.html</filename>
    <includes id="aggregate__across__genes_8hpp" name="aggregate_across_genes.hpp" local="yes" import="no" module="no" objc="no">aggregate_across_genes.hpp</includes>
    <includes id="aggregate__across__cells_8hpp" name="aggregate_across_cells.hpp" local="yes" import="no" module="no" objc="no">aggregate_across_cells.hpp</includes>
    <namespace>scran_aggregate</namespace>
  </compound>
  <compound kind="struct">
    <name>scran_aggregate::AggregateAcrossCellsBuffers</name>
    <filename>structscran__aggregate_1_1AggregateAcrossCellsBuffers.html</filename>
    <templarg>typename Sum_</templarg>
    <templarg>typename Detected_</templarg>
    <templarg>typename Float_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Sum_ * &gt;</type>
      <name>sum</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsBuffers.html</anchorfile>
      <anchor>ab2a2c96b2ee8625591277b4223500d6e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Detected_ * &gt;</type>
      <name>detected</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsBuffers.html</anchorfile>
      <anchor>aecab935899c993b325c6a4666eae882e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Float_ * &gt;</type>
      <name>median</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsBuffers.html</anchorfile>
      <anchor>a25cb5c7e2b874409e3b8c86654e40ac4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_aggregate::AggregateAcrossCellsOptions</name>
    <filename>structscran__aggregate_1_1AggregateAcrossCellsOptions.html</filename>
    <member kind="variable">
      <type>bool</type>
      <name>compute_sum</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsOptions.html</anchorfile>
      <anchor>a3f4836f9e16ff6ffd4c6804541ebc384</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_detected</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsOptions.html</anchorfile>
      <anchor>a698760c10561b3d904cb52abdbbe710e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_median</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsOptions.html</anchorfile>
      <anchor>af4cb0450ae68bbf2287ea0df36413ead</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsOptions.html</anchorfile>
      <anchor>a2469622e7795d0002f6e0b2da4f506f8</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_aggregate::AggregateAcrossCellsResults</name>
    <filename>structscran__aggregate_1_1AggregateAcrossCellsResults.html</filename>
    <templarg>typename Sum_</templarg>
    <templarg>typename Detected_</templarg>
    <templarg>typename Float_</templarg>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Sum_ &gt; &gt;</type>
      <name>sum</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsResults.html</anchorfile>
      <anchor>a06bafa2553d07a80cefd449bea3ab2d9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Detected_ &gt; &gt;</type>
      <name>detected</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsResults.html</anchorfile>
      <anchor>accd409380a04a88b8cf067f1c20889e6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Float_ &gt; &gt;</type>
      <name>median</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsResults.html</anchorfile>
      <anchor>a532ac5295fc04496690c271edfbfed7c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_aggregate::AggregateAcrossGenesBuffers</name>
    <filename>structscran__aggregate_1_1AggregateAcrossGenesBuffers.html</filename>
    <templarg>typename Sum_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Sum_ * &gt;</type>
      <name>sum</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossGenesBuffers.html</anchorfile>
      <anchor>aa739ab1284738c1eaad2c91104c2238e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_aggregate::AggregateAcrossGenesOptions</name>
    <filename>structscran__aggregate_1_1AggregateAcrossGenesOptions.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossGenesOptions.html</anchorfile>
      <anchor>a97c33a8c958774db4b1296b9cd98ddcc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>average</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossGenesOptions.html</anchorfile>
      <anchor>a18c169c320277b91a53d7d4633a3f09e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_aggregate::AggregateAcrossGenesResults</name>
    <filename>structscran__aggregate_1_1AggregateAcrossGenesResults.html</filename>
    <templarg>typename Sum_</templarg>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Sum_ &gt; &gt;</type>
      <name>sum</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossGenesResults.html</anchorfile>
      <anchor>a732d5bac9266655efeaa850f27eb68eb</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_aggregate::AggregateAcrossGenesSet</name>
    <filename>structscran__aggregate_1_1AggregateAcrossGenesSet.html</filename>
    <templarg>typename Gene_</templarg>
    <templarg>typename Weight_</templarg>
    <member kind="function">
      <type></type>
      <name>AggregateAcrossGenesSet</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossGenesSet.html</anchorfile>
      <anchor>a155c8e29c04c17a2a4d25619d1c44aea</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>AggregateAcrossGenesSet</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossGenesSet.html</anchorfile>
      <anchor>afe815c1631fbd39f7b0012eb6e35773f</anchor>
      <arglist>(std::size_t number, const Gene_ *gene, const Weight_ *weight)</arglist>
    </member>
    <member kind="variable">
      <type>std::size_t</type>
      <name>number</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossGenesSet.html</anchorfile>
      <anchor>a7e31865c613fe05f300b14df22fe24ee</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const Gene_ *</type>
      <name>gene</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossGenesSet.html</anchorfile>
      <anchor>a2ebf8fe83566b390f88c90a4c547c9d9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const Weight_ *</type>
      <name>weight</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossGenesSet.html</anchorfile>
      <anchor>a649c6600df0357cf431e36c05f806823</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>scran_aggregate</name>
    <filename>namespacescran__aggregate.html</filename>
    <class kind="struct">scran_aggregate::AggregateAcrossCellsBuffers</class>
    <class kind="struct">scran_aggregate::AggregateAcrossCellsOptions</class>
    <class kind="struct">scran_aggregate::AggregateAcrossCellsResults</class>
    <class kind="struct">scran_aggregate::AggregateAcrossGenesBuffers</class>
    <class kind="struct">scran_aggregate::AggregateAcrossGenesOptions</class>
    <class kind="struct">scran_aggregate::AggregateAcrossGenesResults</class>
    <class kind="struct">scran_aggregate::AggregateAcrossGenesSet</class>
    <member kind="function">
      <type>void</type>
      <name>aggregate_across_cells</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>a4432ed9b60fd5a9e87ea0cc2e3be2825</anchor>
      <arglist>(const tatami::Matrix&lt; Data_, Index_ &gt; &amp;input, const Group_ *const group, const std::size_t num_groups, const AggregateAcrossCellsBuffers&lt; Sum_, Detected_, Float_ &gt; &amp;buffers, const AggregateAcrossCellsOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>AggregateAcrossCellsResults&lt; Sum_, Detected_, Float_ &gt;</type>
      <name>aggregate_across_cells</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>a9602fe5784339e9aaa40aafbc166f775</anchor>
      <arglist>(const tatami::Matrix&lt; Data_, Index_ &gt; &amp;input, const Group_ *const group, const std::size_t num_groups, const AggregateAcrossCellsOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>aggregate_across_genes</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>afaab6fff62f09b24a176f3587ede30fe</anchor>
      <arglist>(const tatami::Matrix&lt; Data_, Index_ &gt; &amp;input, const std::vector&lt; AggregateAcrossGenesSet&lt; Gene_, Weight_ &gt; &gt; &amp;gene_sets, const AggregateAcrossGenesBuffers&lt; Sum_ &gt; &amp;buffers, const AggregateAcrossGenesOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>AggregateAcrossGenesResults&lt; Sum_ &gt;</type>
      <name>aggregate_across_genes</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>a6a74e3bc0960a037637f58a3972cde30</anchor>
      <arglist>(const tatami::Matrix&lt; Data_, Index_ &gt; &amp;input, const std::vector&lt; AggregateAcrossGenesSet&lt; Gene_, Weight_ &gt; &gt; &amp;gene_sets, const AggregateAcrossGenesOptions &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Aggregate expression values across cells</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
