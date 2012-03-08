function barChart() {
  var margin = {top: 0, right: 10, bottom: 20, left: 40},
      width = 960,
      height = 500,
      spacing = .1;

  var key = defaultKey,
      value = defaultValue,
      valueFormat = d3.format(",.f"),
      valueDomain = null;

  var x = d3.scale.linear(),
      y = d3.scale.ordinal(),
      xAxis = d3.svg.axis().scale(x).orient("bottom"),
      yAxis = d3.svg.axis().scale(y).tickSize(0).orient("left");

  function chart(selection) {
    var innerWidth = width - margin.left - margin.right,
        innerHeight = height - margin.top - margin.bottom;

    selection.each(function(data) {

      // Convert the data to standard {key, value} representation.
      // This must be done greedily for nondeterministic accessors.
      data = data.map(function(d, i) {
        return {
          key: "" + key.call(chart, d, i),
          value: +value.call(chart, d, i)
        };
      });

      // Update the x-scale's domain and range.
      // Compute the domain from data if not set explicitly.
      var xExtent = valueDomain || d3.extent(data, defaultValue);
      x   .domain([Math.min(0, xExtent[0]), Math.max(0, xExtent[1])])
          .range([0, innerWidth]);

      // Update the y-scale's domain and range.
      y   .domain(data.map(defaultKey))
          .rangeRoundBands([0, innerHeight], spacing);

      // Stash a snapshot of the new x-scale, and retrieve the old snapshot.
      var x0 = this.__chart__ || x;
      this.__chart__ = x.copy();

      // Select the svg element, if it exists.
      var svg = d3.select(this).selectAll("svg").data([data]);

      // Otherwise, create the skeletal chart.
      var svgEnter = svg.enter().append("svg").append("g");
      svgEnter.append("g").attr("class", "x axis");
      svgEnter.append("g").attr("class", "y axis");

      // Update the outer dimensions.
      var svgUpdate = d3.transition(svg)
          .attr("width", width)
          .attr("height", height);

      // Update the inner dimensions.
      svgUpdate.select("g")
          .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

      // Update the y-axis.
      svgUpdate.select(".x.axis")
          .attr("transform", "translate(0," + innerHeight + ")")
          .call(xAxis);

      // Update the y-axis.
      svgUpdate.select(".y.axis")
          .attr("transform", "translate(" + x(0) + ")")
          .call(yAxis);

      // Update the bars.
      var bar = svg.select("g").selectAll(".bar")
          .data(data, defaultKey);

      // Enter any new bars.
      // Note that we use the new y-scale for the transform;
      // the old y-scale wouldn't have a defined value for the new key.
      var barEnter = bar.enter().insert("g", ".axis")
          .attr("class", "bar")
          .attr("transform", function(d) { return "translate(0," + y(d.key) + ")"; })
          .style("fill-opacity", 1e-6);

      barEnter.append("rect")
          .attr("height", y.rangeBand());

      barEnter.append("text")
          .attr("y", y.rangeBand() / 2)
          .attr("dy", ".35em");

      barEnter.call(position, x0);

      // Enter and update transition.
      var barUpdate = d3.transition(bar)
          .attr("transform", function(d) { return "translate(0," + y(d.key) + ")"; })
          .style("fill-opacity", 1)
          .call(position, x);

      barUpdate.select("rect")
          .attr("height", y.rangeBand());

      barUpdate.select("text")
          .attr("y", y.rangeBand() / 2)
          .text(function(d) { return valueFormat(d.value); });

      // Exit transition.
      var barExit = d3.transition(bar.exit())
          .style("opacity", 1e-6)
          .call(position, x)
          .remove();
    });
  }

  function position(bar, x) {
    bar.select("rect")
        .attr("x", function(d) { return Math.min(x(0), x(d.value)); })
        .attr("width", function(d) { return Math.abs(x(d.value) - x(0)); });

    bar.select("text")
        .attr("x", function(d) { return x(d.value) + (d.value < 0 ? +6 : -6); })
        .attr("text-anchor", function(d) { return d.value < 0 ? "start" : "end"; });
  }

  chart.width = function(_) {
    if (!arguments.length) return width;
    width = _;
    return chart;
  };

  chart.height = function(_) {
    if (!arguments.length) return height;
    height = _;
    return chart;
  };

  chart.margin = function(_) {
    if (!arguments.length) return margin;
    margin = _;
    return chart;
  };

  chart.spacing = function(_) {
    if (!arguments.length) return spacing;
    spacing = _;
    return chart;
  };

  chart.key = function(_) {
    if (!arguments.length) return key;
    key = _;
    return chart;
  };

  chart.keyFormat = function(_) {
    if (!arguments.length) return yAxis.tickFormat();
    yAxis.tickFormat(_);
    return chart;
  };

  chart.value = function(_) {
    if (!arguments.length) return value;
    value = _;
    return chart;
  };

  chart.valueFormat = function(_) {
    if (!arguments.length) return valueFormat;
    valueFormat = _;
    return chart;
  };

  chart.valueDomain = function(_) {
    if (!arguments.length) return valueDomain;
    valueDomain = _;
    return chart;
  };

  function defaultKey(d) { return d.key; }
  function defaultValue(d) { return d.value; }

  return d3.rebind(chart, xAxis, "ticks", "tickFormat", "tickSize", "tickPadding");
}
