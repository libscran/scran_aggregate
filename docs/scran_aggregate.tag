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
    <class kind="struct">scran_aggregate::AggregateAcrossGenesResults</class>
    <namespace>scran_aggregate</namespace>
  </compound>
  <compound kind="file">
    <name>clean_factor.hpp</name>
    <path>scran_aggregate/</path>
    <filename>clean__factor_8hpp.html</filename>
    <namespace>scran_aggregate</namespace>
  </compound>
  <compound kind="file">
    <name>combine_factors.hpp</name>
    <path>scran_aggregate/</path>
    <filename>combine__factors_8hpp.html</filename>
    <includes id="clean__factor_8hpp" name="clean_factor.hpp" local="yes" import="no" module="no" objc="no">clean_factor.hpp</includes>
    <namespace>scran_aggregate</namespace>
  </compound>
  <compound kind="file">
    <name>scran_aggregate.hpp</name>
    <path>scran_aggregate/</path>
    <filename>scran__aggregate_8hpp.html</filename>
    <includes id="aggregate__across__genes_8hpp" name="aggregate_across_genes.hpp" local="yes" import="no" module="no" objc="no">aggregate_across_genes.hpp</includes>
    <includes id="aggregate__across__cells_8hpp" name="aggregate_across_cells.hpp" local="yes" import="no" module="no" objc="no">aggregate_across_cells.hpp</includes>
    <includes id="combine__factors_8hpp" name="combine_factors.hpp" local="yes" import="no" module="no" objc="no">combine_factors.hpp</includes>
    <includes id="clean__factor_8hpp" name="clean_factor.hpp" local="yes" import="no" module="no" objc="no">clean_factor.hpp</includes>
    <namespace>scran_aggregate</namespace>
  </compound>
  <compound kind="struct">
    <name>scran_aggregate::AggregateAcrossCellsBuffers</name>
    <filename>structscran__aggregate_1_1AggregateAcrossCellsBuffers.html</filename>
    <templarg>typename Sum_</templarg>
    <templarg>typename Detected_</templarg>
    <member kind="variable">
      <type>std::vector&lt; Sum_ * &gt;</type>
      <name>sums</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsBuffers.html</anchorfile>
      <anchor>af18936a8bb6ba1fbdb32f4584683288c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; Detected_ * &gt;</type>
      <name>detected</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsBuffers.html</anchorfile>
      <anchor>a468273e5a9566bbc2caae8e204a0d070</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran_aggregate::AggregateAcrossCellsOptions</name>
    <filename>structscran__aggregate_1_1AggregateAcrossCellsOptions.html</filename>
    <member kind="variable">
      <type>bool</type>
      <name>compute_sums</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsOptions.html</anchorfile>
      <anchor>ad45ae84870be13fad30a9dbdd0299d62</anchor>
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
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Sum_ &gt; &gt;</type>
      <name>sums</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsResults.html</anchorfile>
      <anchor>ad09fedfa434a9be5074b47fd950098b7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Detected_ &gt; &gt;</type>
      <name>detected</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsResults.html</anchorfile>
      <anchor>abd1392875a4be38f5a7317b128c7b13b</anchor>
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
  <compound kind="namespace">
    <name>scran_aggregate</name>
    <filename>namespacescran__aggregate.html</filename>
    <class kind="struct">scran_aggregate::AggregateAcrossCellsBuffers</class>
    <class kind="struct">scran_aggregate::AggregateAcrossCellsOptions</class>
    <class kind="struct">scran_aggregate::AggregateAcrossCellsResults</class>
    <class kind="struct">scran_aggregate::AggregateAcrossGenesBuffers</class>
    <class kind="struct">scran_aggregate::AggregateAcrossGenesOptions</class>
    <class kind="struct">scran_aggregate::AggregateAcrossGenesResults</class>
    <member kind="function">
      <type>void</type>
      <name>aggregate_across_cells</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>a455854154026f84c3d8c52a216472c2b</anchor>
      <arglist>(const tatami::Matrix&lt; Data_, Index_ &gt; &amp;input, const Group_ *const group, const AggregateAcrossCellsBuffers&lt; Sum_, Detected_ &gt; &amp;buffers, const AggregateAcrossCellsOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>AggregateAcrossCellsResults&lt; Sum_, Detected_ &gt;</type>
      <name>aggregate_across_cells</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>a14dbda14b9307f2c4dbedfd8c3622f5c</anchor>
      <arglist>(const tatami::Matrix&lt; Data_, Index_ &gt; &amp;input, const Group_ *const group, const AggregateAcrossCellsOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>aggregate_across_genes</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>abf6732ea8012e17cc65039228663eae0</anchor>
      <arglist>(const tatami::Matrix&lt; Data_, Index_ &gt; &amp;input, const std::vector&lt; std::tuple&lt; std::size_t, const Gene_ *, const Weight_ * &gt; &gt; &amp;gene_sets, const AggregateAcrossGenesBuffers&lt; Sum_ &gt; &amp;buffers, const AggregateAcrossGenesOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>AggregateAcrossGenesResults&lt; Sum_ &gt;</type>
      <name>aggregate_across_genes</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>a20815db8de3133fc282433dfa74443d9</anchor>
      <arglist>(const tatami::Matrix&lt; Data_, Index_ &gt; &amp;input, const std::vector&lt; std::tuple&lt; std::size_t, const Gene_ *, const Weight_ * &gt; &gt; &amp;gene_sets, const AggregateAcrossGenesOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; Factor_ &gt;</type>
      <name>clean_factor</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>ac14284076a7d9332cd95646e8db4e769</anchor>
      <arglist>(const std::size_t n, const Factor_ *const factor, Output_ *const cleaned)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::vector&lt; Factor_ &gt; &gt;</type>
      <name>combine_factors</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>ad80f3b36beb047078b69038867e8e714</anchor>
      <arglist>(const std::size_t n, const std::vector&lt; const Factor_ * &gt; &amp;factors, Combined_ *const combined)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; std::vector&lt; Factor_ &gt; &gt;</type>
      <name>combine_factors_unused</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>a252b9226f106b365c068f50a206597e5</anchor>
      <arglist>(const std::size_t n, const std::vector&lt; std::pair&lt; const Factor_ *, Number_ &gt; &gt; &amp;factors, Combined_ *const combined)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Aggregate expression values across cells</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
