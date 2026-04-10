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
      <name>sums</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsBuffers.html</anchorfile>
      <anchor>ac4e883670838418b0a61bb7f254b7065</anchor>
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
      <name>medians</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsBuffers.html</anchorfile>
      <anchor>a23d5b96d6ca63bf1fcc6d5c9ea823226</anchor>
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
      <type>bool</type>
      <name>compute_medians</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsOptions.html</anchorfile>
      <anchor>a4c52daf55fc6064926d8ae5eae624480</anchor>
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
      <name>sums</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsResults.html</anchorfile>
      <anchor>aaffd51d13ca7681a494b904d09019ef4</anchor>
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
      <name>medians</name>
      <anchorfile>structscran__aggregate_1_1AggregateAcrossCellsResults.html</anchorfile>
      <anchor>aa2ac90592c80f617ce2052d2684e5613</anchor>
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
      <anchor>ace4c8f8b87106e77dedd5d394fa370f1</anchor>
      <arglist>(const tatami::Matrix&lt; Data_, Index_ &gt; &amp;input, const Group_ *const group, const AggregateAcrossCellsBuffers&lt; Sum_, Detected_, Float_ &gt; &amp;buffers, const AggregateAcrossCellsOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>AggregateAcrossCellsResults&lt; Sum_, Detected_, Float_ &gt;</type>
      <name>aggregate_across_cells</name>
      <anchorfile>namespacescran__aggregate.html</anchorfile>
      <anchor>a56acd20f25216e5149635077ea76732b</anchor>
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
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Aggregate expression values across cells</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
