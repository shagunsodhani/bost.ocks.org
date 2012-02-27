function emptyChart() {
  var width = 960,
      height = 120;

  function chart(div) {
    div.each(function() {
      var div = d3.select(this);

      // Select the vg element, if it exists.
      var svg = div.selectAll("svg").data([null]);

      // Otherwise, create the svg element.
      svg.enter().append("svg");

      // Update the outer dimensions.
      svg .attr("width", width)
          .attr("height", height);
    });
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

  return chart;
}

function marginChart() {
  var margin = {top: 0, right: 0, bottom: 0, left: 0},
      width = 960,
      height = 120;

  function chart(div) {
    div.each(function() {
      var div = d3.select(this);

      // Select the svg element, if it exists.
      var svg = div.selectAll("svg").data([null]);

      // Otherwise, create the svg and rect elements.
      svg.enter().append("svg").append("rect");

      // Update the outer dimensions.
      svg .attr("width", width)
          .attr("height", height);

      // Update the inner dimensions.
      svg.select("rect")
          .attr("x", margin.left)
          .attr("y", margin.top)
          .attr("width", width - margin.left - margin.right)
          .attr("height", height - margin.top - margin.bottom);
    });
  }

  chart.margin = function(_) {
    if (!arguments.length) return margin;
    margin = _;
    return chart;
  };

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

  return chart;
}

function axisChart() {
  var margin = {top: 0, right: 0, bottom: 20, left: 0},
      width = 960,
      height = 120,
      x = d3.scale.linear(),
      xAxis = d3.svg.axis().scale(x).orient("bottom");

  function chart(div) {
    div.each(function() {
      var div = d3.select(this);

      // Select the svg element, if it exists.
      var svg = div.selectAll("svg").data([null]);

      // Otherwise, create the skeletal chart.
      var gEnter = svg.enter().append("svg").append("g");
      gEnter.append("g").attr("class", "x axis");

      // Update the outer dimensions.
      svg .attr("width", width)
          .attr("height", height);

      // Update the inner dimensions.
      var g = svg.select("g")
          .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

      // Compute the inner dimensions.
      var innerWidth = width - margin.left - margin.right,
          innerHeight = height - margin.top - margin.bottom;

      // Update the x-scale.
      x.range([0, innerWidth]);

      // Update the x-axis.
      g.select(".x.axis")
          .attr("transform", "translate(0," + innerHeight + ")")
          .call(xAxis);
    });
  }

  chart.margin = function(_) {
    if (!arguments.length) return margin;
    margin = _;
    return chart;
  };

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

  return chart;
}

function areaChart() {
  var margin = {top: 0, right: 0, bottom: 20, left: 0},
      width = 960,
      height = 120,
      x = d3.scale.linear(),
      y = d3.scale.linear(),
      xAxis = d3.svg.axis().scale(x).orient("bottom"),
      area = d3.svg.area().x(X).y1(Y);

  function chart(div) {
    div.each(function(d) {
      var div = d3.select(this);

      // Select the svg element, if it exists.
      var svg = div.selectAll("svg").data([null]);

      // Otherwise, create the skeletal chart.
      var gEnter = svg.enter().append("svg").append("g");
      gEnter.append("path").attr("class", "area");
      gEnter.append("g").attr("class", "x axis");

      // Update the outer dimensions.
      svg .attr("width", width)
          .attr("height", height);

      // Update the inner dimensions.
      var g = svg.select("g")
          .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

      // Compute the inner dimensions.
      x.range([0, width - margin.left - margin.right]);
      y.range([height - margin.top - margin.bottom, 0]);

      // Update the path.
      g.select(".area").attr("d", area.y0(y.range()[0])(d));

      // Update the x-axis.
      g.select(".x.axis")
          .attr("transform", "translate(0," + y.range()[0] + ")")
          .call(xAxis);
    });
  }

  function X(d) {
    return x(d[0]);
  }

  function Y(d) {
    return y(d[1]);
  }

  chart.margin = function(_) {
    if (!arguments.length) return margin;
    margin = _;
    return chart;
  };

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

  return chart;
}
