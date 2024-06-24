<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.8">
  <compound kind="file">
    <name>aggregate_across_cells.hpp</name>
    <path>scran/</path>
    <filename>aggregate__across__cells_8hpp.html</filename>
    <class kind="struct">scran::aggregate_across_cells::Options</class>
    <class kind="struct">scran::aggregate_across_cells::Results</class>
    <namespace>scran</namespace>
    <namespace>scran::aggregate_across_cells</namespace>
  </compound>
  <compound kind="file">
    <name>combine_factors.hpp</name>
    <path>scran/</path>
    <filename>combine__factors_8hpp.html</filename>
    <class kind="struct">scran::combine_factors::Results</class>
    <namespace>scran</namespace>
    <namespace>scran::combine_factors</namespace>
  </compound>
  <compound kind="file">
    <name>scran.hpp</name>
    <path>scran/</path>
    <filename>scran_8hpp.html</filename>
    <namespace>scran</namespace>
  </compound>
  <compound kind="struct">
    <name>scran::aggregate_across_cells::Options</name>
    <filename>structscran_1_1aggregate__across__cells_1_1Options.html</filename>
    <member kind="variable">
      <type>bool</type>
      <name>compute_sums</name>
      <anchorfile>structscran_1_1aggregate__across__cells_1_1Options.html</anchorfile>
      <anchor>abafe797e3f1808a7d258aeed7e1fe823</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>compute_detected</name>
      <anchorfile>structscran_1_1aggregate__across__cells_1_1Options.html</anchorfile>
      <anchor>a27912761546ec8d507041c15b648b561</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran_1_1aggregate__across__cells_1_1Options.html</anchorfile>
      <anchor>a7cec296b4ee8633c76d40ca0d6c23666</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran::aggregate_across_cells::Results</name>
    <filename>structscran_1_1aggregate__across__cells_1_1Results.html</filename>
    <templarg>typename Sum_</templarg>
    <templarg>typename Detected_</templarg>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Sum_ &gt; &gt;</type>
      <name>sums</name>
      <anchorfile>structscran_1_1aggregate__across__cells_1_1Results.html</anchorfile>
      <anchor>a7df8887bf412f7f3fb0fd08becb42da1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Detected_ &gt; &gt;</type>
      <name>detected</name>
      <anchorfile>structscran_1_1aggregate__across__cells_1_1Results.html</anchorfile>
      <anchor>a5f490db8a7f5e5408e14235f28b775c6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran::combine_factors::Results</name>
    <filename>structscran_1_1combine__factors_1_1Results.html</filename>
    <templarg>typename Factor_</templarg>
    <member kind="variable">
      <type>std::vector&lt; std::vector&lt; Factor_ &gt; &gt;</type>
      <name>factors</name>
      <anchorfile>structscran_1_1combine__factors_1_1Results.html</anchorfile>
      <anchor>a2ab4dea8b028bb3cf539c9f14f799f7c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; size_t &gt;</type>
      <name>counts</name>
      <anchorfile>structscran_1_1combine__factors_1_1Results.html</anchorfile>
      <anchor>a0eb2102be859441545431e2470c3109d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>scran</name>
    <filename>namespacescran.html</filename>
    <namespace>scran::aggregate_across_cells</namespace>
    <namespace>scran::combine_factors</namespace>
  </compound>
  <compound kind="namespace">
    <name>scran::aggregate_across_cells</name>
    <filename>namespacescran_1_1aggregate__across__cells.html</filename>
    <class kind="struct">scran::aggregate_across_cells::Options</class>
    <class kind="struct">scran::aggregate_across_cells::Results</class>
    <member kind="function">
      <type>void</type>
      <name>compute</name>
      <anchorfile>namespacescran_1_1aggregate__across__cells.html</anchorfile>
      <anchor>aa8c09dfee9f424daa59e29f2826460dc</anchor>
      <arglist>(const tatami::Matrix&lt; Data_, Index_ &gt; *input, const Factor_ *factor, std::vector&lt; Sum_ * &gt; sums, std::vector&lt; Detected_ * &gt; detected, const Options &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>Results&lt; Sum_, Detected_ &gt;</type>
      <name>compute</name>
      <anchorfile>namespacescran_1_1aggregate__across__cells.html</anchorfile>
      <anchor>acd5ae1de15376ecfeae0baf4030366a8</anchor>
      <arglist>(const tatami::Matrix&lt; Data_, Index_ &gt; *input, const Factor_ *factor, const Options &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>scran::combine_factors</name>
    <filename>namespacescran_1_1combine__factors.html</filename>
    <class kind="struct">scran::combine_factors::Results</class>
    <member kind="function">
      <type>Results&lt; Factor_ &gt;</type>
      <name>compute</name>
      <anchorfile>namespacescran_1_1combine__factors.html</anchorfile>
      <anchor>a522c70000c2235e998e112489b849f02</anchor>
      <arglist>(size_t n, const std::vector&lt; const Factor_ * &gt; &amp;factors, Combined_ *combined)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Aggregate expression values across cells</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
