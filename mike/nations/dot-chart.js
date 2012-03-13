function dotChart() {
  var margin = {top: 19.5, right: 19.5, bottom: 19.5, left: 39.5},
      width = 960,
      height = 500;

  var xValue = function(d) { return d[0]; },
      yValue = function(d) { return d[1]; },
      xDomain = null,
      yDomain = null,
      xAxis = d3.svg.axis().orient("bottom"),
      yAxis = d3.svg.axis().orient("left");

  var rValue = function() { return 1; },
      r = d3.scale.sqrt();

  var zValue = function() { return "undefined"; },
      z = d3.scale.category10();

  var keyValue = function(d, i) { return i; };

  function chart(selection) {
    var innerWidth = width - margin.left - margin.right,
        innerHeight = height - margin.top - margin.bottom,
        x = xAxis.scale(),
        y = yAxis.scale();

    // Update the scales' ranges.
    x.range([0, innerWidth]);
    y.range([innerHeight, 0]);
    r.range([0, Math.min(innerWidth, innerHeight) / 12]);

    selection.each(function(data) {

      // Convert the data to standard representation.
      // This must be done greedily for nondeterministic accessors.
      data = data.map(function(d, i) {
        return {
          data: d,
          index: i,
          x: +xValue.call(chart, d, i),
          y: +yValue.call(chart, d, i),
          radius: +rValue.call(chart, d, i),
          color: "" + zValue.call(chart, d, i),
          key: "" + keyValue.call(chart, d, i)
        };
      });

      // Update the scales' domains.
      // TODO Allow the color scale to be quantitative (use extent, not map).
      x   .domain(xDomain || d3.extent(data, function(d) { return d.x; }));
      y   .domain(yDomain || d3.extent(data, function(d) { return d.y; }));
      r   .domain([0, d3.max(data, function(d) { return d.radius; })]);
      z   .domain(data.map(function(d) { return d.color; }));

      // Stash a snapshot of the new scales, and retrieve the old snapshot.
      var snapshot = this.__chart__, x0 = x, y0 = y, r0 = r;
      if (snapshot) x0 = snapshot.x, y0 = snapshot.y, r0 = snapshot.radius;
      this.__chart__ = {x: x.copy(), y: y.copy(), radius: r.copy()};

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

      // Update the dots.
      var dot = svg.select("g").selectAll(".dot")
          .data(data, function(d) { return d.key; });

      // Enter any new dots.
      var dotEnter = dot.enter().append("circle")
          .attr("class", "dot")
          .style("stroke-opacity", 1e-6)
          .style("fill-opacity", 1e-6)
          .call(encode, x0, y0, r0);

      // Sort the dots by descending radius.
      dot.sort(function(a, b) { return b.radius - a.radius; });

      // Enter and update transition.
      var dotUpdate = d3.transition(dot)
          .style("stroke-opacity", 1)
          .style("fill-opacity", 1)
          .call(encode, x, y, r);

      // Exit transition.
      var dotExit = d3.transition(dot.exit())
          .style("stroke-opacity", 1e-6)
          .style("fill-opacity", 1e-6)
          .call(encode, x, y, r)
          .remove();
    });
  }

  function encode(dot, x, y, r) {
    dot .attr("cx", function(d) { return x(d.x); })
        .attr("cy", function(d) { return y(d.y); })
        .attr("r", function(d) { return r(d.radius); })
        .style("fill", function(d) { return z(d.color); });
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

  chart.key = function(_) {
    if (!arguments.length) return keyValue;
    keyValue = _;
    return chart;
  };

  chart.x = function(_) {
    if (!arguments.length) return xValue;
    xValue = _;
    return chart;
  };

  chart.xDomain = function(_) {
    if (!arguments.length) return xDomain;
    xDomain = _;
    return chart;
  };

  chart.xAxis = function(_) {
    if (!arguments.length) return xAxis;
    xAxis = _;
    return chart;
  };

  chart.y = function(_) {
    if (!arguments.length) return yValue;
    yValue = _;
    return chart;
  };

  chart.yDomain = function(_) {
    if (!arguments.length) return yDomain;
    yDomain = _;
    return chart;
  };

  chart.yAxis = function(_) {
    if (!arguments.length) return yAxis;
    yAxis = _;
    return chart;
  };

  chart.radius = function(_) {
    if (!arguments.length) return rValue;
    rValue = _;
    return chart;
  };

  chart.color = function(_) {
    if (!arguments.length) return zValue;
    zValue = _;
    return chart;
  };

  chart.colorScale = function(_) {
    if (!arguments.length) return z;
    z = _;
    return chart;
  };

  return chart;
}
