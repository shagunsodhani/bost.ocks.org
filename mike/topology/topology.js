function topology(objects, Q, callback) {
  if (arguments.length < 3) callback = Q, Q = 0;

  var hashBuffer = new ArrayBuffer(8),
      hashFloats = new Float64Array(hashBuffer),
      hashInts = new Int32Array(hashBuffer);

  function hashFloat(x) {
    hashFloats[0] = x;
    x = hashInts[1] ^ hashInts[0];
    x ^= (x >>> 20) ^ (x >>> 12);
    x ^= (x >>> 7) ^ (x >>> 4);
    return x;
  }

  function hashPoint(point) {
    var h = (hashFloat(point[0]) + 31 * hashFloat(point[1])) | 0;
    return h < 0 ? ~h : h;
  }

  function equalPoint(pointA, pointB) {
    return pointA[0] === pointB[0] && pointA[1] === pointB[1];
  }

  function bounds(objects) {
    var x0 = Infinity,
        y0 = Infinity,
        x1 = -Infinity,
        y1 = -Infinity;

    function boundGeometry(geometry) {
      if (geometry && boundGeometryType.hasOwnProperty(geometry.type)) boundGeometryType[geometry.type](geometry);
    }

    var boundGeometryType = {
      GeometryCollection: function(o) { o.geometries.forEach(boundGeometry); },
      Point: function(o) { boundPoint(o.coordinates); },
      MultiPoint: function(o) { o.coordinates.forEach(boundPoint); },
      LineString: function(o) { boundLine(o.coordinates); },
      MultiLineString: function(o) { o.coordinates.forEach(boundLine); },
      Polygon: function(o) { o.coordinates.forEach(boundLine); },
      MultiPolygon: function(o) { o.coordinates.forEach(boundMultiLine); }
    };

    function boundPoint(coordinates) {
      var x = coordinates[0],
          y = coordinates[1];
      if (x < x0) x0 = x;
      if (x > x1) x1 = x;
      if (y < y0) y0 = y;
      if (y > y1) y1 = y;
    }

    function boundLine(coordinates) {
      coordinates.forEach(boundPoint);
    }

    function boundMultiLine(coordinates) {
      coordinates.forEach(boundLine);
    }

    for (var key in objects) {
      boundGeometry(objects[key]);
    }

    return [x0, y0, x1, y1];
  }

  function quantize(objects, bbox, Q) {
    var x0 = isFinite(bbox[0]) ? bbox[0] : 0,
        y0 = isFinite(bbox[1]) ? bbox[1] : 0,
        x1 = isFinite(bbox[2]) ? bbox[2] : 0,
        y1 = isFinite(bbox[3]) ? bbox[3] : 0,
        kx = x1 - x0 ? (Q - 1) / (x1 - x0) : 1,
        ky = y1 - y0 ? (Q - 1) / (y1 - y0) : 1;

    function quantizeGeometry(geometry) {
      if (geometry && quantizeGeometryType.hasOwnProperty(geometry.type)) quantizeGeometryType[geometry.type](geometry);
    }

    var quantizeGeometryType = {
      GeometryCollection: function(o) { o.geometries.forEach(quantizeGeometry); },
      Point: function(o) { quantizePoint(o.coordinates); },
      MultiPoint: function(o) { o.coordinates.forEach(quantizePoint); },
      LineString: function(o) {
        var line = o.coordinates;
        quantizeLine(line);
        if (line.length < 2) line[1] = line[0]; // must have 2+
      },
      MultiLineString: function(o) {
        for (var lines = o.coordinates, i = 0, n = lines.length; i < n; ++i) {
          var line = lines[i];
          quantizeLine(line);
          if (line.length < 2) line[1] = line[0]; // must have 2+
        }
      },
      Polygon: function(o) {
        for (var rings = o.coordinates, i = 0, n = rings.length; i < n; ++i) {
          var ring = rings[i];
          quantizeLine(ring);
          while (ring.length < 4) ring.push(ring[0]); // must have 4+
        }
      },
      MultiPolygon: function(o) {
        for (var polygons = o.coordinates, i = 0, n = polygons.length; i < n; ++i) {
          for (var rings = polygons[i], j = 0, m = rings.length; j < m; ++j) {
            var ring = rings[j];
            quantizeLine(ring);
            while (ring.length < 4) ring.push(ring[0]); // must have 4+
          }
        }
      }
    };

    function quantizePoint(coordinates) {
      coordinates[0] = Math.round((coordinates[0] - x0) * kx);
      coordinates[1] = Math.round((coordinates[1] - y0) * ky);
    }

    function quantizeLine(coordinates) {
      var i = 0,
          j = 1,
          n = coordinates.length,
          pi = coordinates[0],
          pj,
          px = pi[0] = Math.round((pi[0] - x0) * kx),
          py = pi[1] = Math.round((pi[1] - y0) * ky),
          x,
          y;

      while (++i < n) {
        pi = coordinates[i];
        x = Math.round((pi[0] - x0) * kx);
        y = Math.round((pi[1] - y0) * ky);
        if (x !== px || y !== py) { // skip coincident points
          pj = coordinates[j++];
          pj[0] = px = x;
          pj[1] = py = y;
        }
      }

      coordinates.length = j;
    }

    for (var key in objects) {
      quantizeGeometry(objects[key]);
    }

    return {
      scale: [1 / kx, 1 / ky],
      translate: [x0, y0]
    };
  }


  function extract(objects) {
    var index = -1,
        lines = [],
        rings = [],
        coordinates = [];

    function extractGeometry(geometry) {
      if (geometry && extractGeometryType.hasOwnProperty(geometry.type)) extractGeometryType[geometry.type](geometry);
    }

    var extractGeometryType = {
      GeometryCollection: function(o) { o.geometries.forEach(extractGeometry); },
      LineString: function(o) { o.arcs = extractLine(o.coordinates); delete o.coordinates; },
      MultiLineString: function(o) { o.arcs = o.coordinates.map(extractLine); delete o.coordinates; },
      Polygon: function(o) { o.arcs = o.coordinates.map(extractRing); delete o.coordinates; },
      MultiPolygon: function(o) { o.arcs = o.coordinates.map(extractMultiRing); delete o.coordinates; }
    };

    function extractLine(line) {
      for (var i = 0, n = line.length; i < n; ++i) coordinates[++index] = line[i];
      var arc = {0: index - n + 1, 1: index};
      lines.push(arc);
      callback({type: "extract:line", coordinates: coordinates.slice(index - n + 1, index + 1)});
      return arc;
    }

    function extractRing(ring) {
      for (var i = 0, n = ring.length; i < n; ++i) coordinates[++index] = ring[i];
      var arc = {0: index - n + 1, 1: index};
      rings.push(arc);
      callback({type: "extract:ring", coordinates: coordinates.slice(index - n + 1, index + 1)});
      return arc;
    }

    function extractMultiRing(rings) {
      return rings.map(extractRing);
    }

    for (var key in objects) {
      extractGeometry(objects[key]);
    }

    return {
      type: "Topology",
      coordinates: coordinates,
      lines: lines,
      rings: rings,
      objects: objects
    };
  }

  function join(topology) {
    var coordinates = topology.coordinates,
        lines = topology.lines,
        rings = topology.rings,
        visitedByPoint,
        neighborsByPoint = hashtable(coordinates.length, hashPoint, equalPoint),
        junctionByPoint = hashtable(coordinates.length, hashPoint, equalPoint);

    for (var i = 0, n = lines.length; i < n; ++i) {
      var line = lines[i],
          lineStart = line[0],
          lineEnd = line[1],
          previousPoint = null,
          currentPoint = coordinates[lineStart],
          nextPoint = coordinates[++lineStart];
      visitedByPoint = hashtable(lineEnd - lineStart, hashPoint, equalPoint);
      junctionByPoint.set(currentPoint, true), callback({type: "join:end", coordinates: currentPoint}); // start
      while (++lineStart <= lineEnd) {
        sequence(previousPoint = currentPoint, currentPoint = nextPoint, nextPoint = coordinates[lineStart]);
      }
      junctionByPoint.set(nextPoint, true), callback({type: "join:end", coordinates: nextPoint}); // end
    }

    for (var i = 0, n = rings.length; i < n; ++i) {
      var ring = rings[i],
          ringStart = ring[0] + 1,
          ringEnd = ring[1],
          previousPoint = coordinates[ringEnd - 1],
          currentPoint = coordinates[ringStart - 1],
          nextPoint = coordinates[ringStart];
      visitedByPoint = hashtable(ringEnd - ringStart + 1, hashPoint, equalPoint);
      sequence(previousPoint, currentPoint, nextPoint);
      while (++ringStart <= ringEnd) {
        sequence(previousPoint = currentPoint, currentPoint = nextPoint, nextPoint = coordinates[ringStart]);
      }
    }

    function sequence(previousPoint, currentPoint, nextPoint) {
      if (visitedByPoint.get(currentPoint)) return;
      visitedByPoint.set(currentPoint, true);
      var neighbors = neighborsByPoint.get(currentPoint);
      if (neighbors) {
        if (!(equalPoint(neighbors[0], previousPoint)
          && equalPoint(neighbors[1], nextPoint))
          && !(equalPoint(neighbors[0], nextPoint)
          && equalPoint(neighbors[1], previousPoint))) {
          if (!junctionByPoint.get(currentPoint)) junctionByPoint.set(currentPoint, true), callback({type: "join:middle", coordinates: [[previousPoint, currentPoint, nextPoint], [neighbors[0], currentPoint, neighbors[1]]]});
        }
      } else {
        neighborsByPoint.set(currentPoint, [previousPoint, nextPoint]);
      }
    }

    return junctionByPoint;
  };

  function hashtable(size, hash, equal) {
    var hashtable = new Array(size = 1 << Math.ceil(Math.log(size + 1) / Math.LN2)),
        mask = size - 1;

    function set(key, value) {
      var index = hash(key) & mask,
          match = hashtable[index];
      while (match != null) {
        if (equal(match.key, key)) return match.value = value;
        match = hashtable[index = (index + 1) & mask];
      }
      hashtable[index] = {key: key, value: value};
      return value;
    }

    function get(key, missingValue) {
      var index = hash(key) & mask,
          match = hashtable[index];
      while (match != null) {
        if (equal(match.key, key)) return match.value;
        match = hashtable[index = (index + 1) & mask];
      }
      return missingValue;
    }

    function remove(key) {
      var index = hash(key) & mask,
          match = hashtable[index];
      while (match != null) {
        if (equal(match.key, key)) {
          hashtable[index] = null;
          match = hashtable[index = (index + 1) & mask];
          if (match != null) {
            hashtable[index] = null;
            set(match.key, match.value);
          }
          return true;
        }
        match = hashtable[index = (index + 1) & mask];
      }
      return false;
    }

    function keys() {
      var keys = [];
      for (var i = 0, n = hashtable.length; i < n; ++i) {
        var match = hashtable[i];
        if (match != null) keys.push(match.key);
      }
      return keys;
    }

    return {
      set: set,
      get: get,
      remove: remove,
      keys: keys
    };
  }

  function dedup(topology) {
    var coordinates = topology.coordinates,
        lines = topology.lines,
        rings = topology.rings,
        arcCount = lines.length + rings.length;

    delete topology.lines;
    delete topology.rings;

    // Count the number of (non-unique) arcs to initialize the hashtable safely.
    for (var i = 0, n = lines.length; i < n; ++i) {
      var line = lines[i]; while (line = line.next) ++arcCount;
    }
    for (var i = 0, n = rings.length; i < n; ++i) {
      var ring = rings[i]; while (ring = ring.next) ++arcCount;
    }

    var arcsByEnd = hashtable(arcCount * 2, hashPoint, equalPoint),
        arcs = topology.arcs = [];

    for (var i = 0, n = lines.length; i < n; ++i) {
      var line = lines[i];
      do {
        dedupLine(line);
      } while (line = line.next);
    }

    for (var i = 0, n = rings.length; i < n; ++i) {
      var ring = rings[i];
      if (ring.next) {
        do {
          dedupLine(ring);
        } while (ring = ring.next);
      } else {
        dedupRing(ring);
      }
    }

    function dedupLine(arc) {
      var startPoint,
          endPoint,
          startArcs,
          endArcs;

      if (startArcs = arcsByEnd.get(startPoint = coordinates[arc[0]])) {
        for (var i = 0, n = startArcs.length; i < n; ++i) {
          var startArc = startArcs[i];
          if (equalLine(startArc, arc)) {
            callback({type: "dedup:oldline", coordinates: coordinates.slice(startArc[0], startArc[1] + 1)});
            arc[0] = startArc[0];
            arc[1] = startArc[1];
            return;
          }
        }
      }

      if (endArcs = arcsByEnd.get(endPoint = coordinates[arc[1]])) {
        for (var i = 0, n = endArcs.length; i < n; ++i) {
          var endArc = endArcs[i];
          if (reverseEqualLine(endArc, arc)) {
            callback({type: "dedup:oldline", coordinates: coordinates.slice(endArc[0], endArc[1] + 1)});
            arc[1] = endArc[0];
            arc[0] = endArc[1];
            return;
          }
        }
      }

      callback({type: "dedup:newline", coordinates: coordinates.slice(arc[0], arc[1] + 1)});
      if (startArcs) startArcs.push(arc); else arcsByEnd.set(startPoint, [arc]);
      if (endArcs) endArcs.push(arc); else arcsByEnd.set(endPoint, [arc]);
      arcs.push(arc);
    }

    function dedupRing(arc) {
      var endPoint,
          endArcs;

      if (endArcs = arcsByEnd.get(endPoint = coordinates[arc[0]])) {
        for (var i = 0, n = endArcs.length; i < n; ++i) {
          var endArc = endArcs[i];
          if (equalRing(endArc, arc)) {
            callback({type: "dedup:oldring", coordinates: coordinates.slice(endArc[0], endArc[1] + 1)});
            arc[0] = endArc[0];
            arc[1] = endArc[1];
            return;
          }
          if (reverseEqualRing(endArc, arc)) {
            callback({type: "dedup:oldring", coordinates: coordinates.slice(endArc[0], endArc[1] + 1)});
            arc[0] = endArc[1];
            arc[1] = endArc[0];
            return;
          }
        }
      }

      if (endArcs = arcsByEnd.get(endPoint = coordinates[arc[0] + findMinimumOffset(arc)])) {
        for (var i = 0, n = endArcs.length; i < n; ++i) {
          var endArc = endArcs[i];
          if (equalRing(endArc, arc)) {
            callback({type: "dedup:oldring", coordinates: coordinates.slice(endArc[0], endArc[1] + 1)});
            arc[0] = endArc[0];
            arc[1] = endArc[1];
            return;
          }
          if (reverseEqualRing(endArc, arc)) {
            callback({type: "dedup:oldring", coordinates: coordinates.slice(endArc[0], endArc[1] + 1)});
            arc[0] = endArc[1];
            arc[1] = endArc[0];
            return;
          }
        }
      }

      callback({type: "dedup:newring", coordinates: coordinates.slice(arc[0], arc[1] + 1)});
      if (endArcs) endArcs.push(arc); else arcsByEnd.set(endPoint, [arc]);
      arcs.push(arc);
    }

    function equalLine(arcA, arcB) {
      var ia = arcA[0], ib = arcB[0],
          ja = arcA[1], jb = arcB[1];
      if (ia - ja !== ib - jb) return false;
      for (; ia <= ja; ++ia, ++ib) if (!equalPoint(coordinates[ia], coordinates[ib])) return false;
      return true;
    }

    function reverseEqualLine(arcA, arcB) {
      var ia = arcA[0], ib = arcB[0],
          ja = arcA[1], jb = arcB[1];
      if (ia - ja !== ib - jb) return false;
      for (; ia <= ja; ++ia, --jb) if (!equalPoint(coordinates[ia], coordinates[jb])) return false;
      return true;
    }

    function equalRing(arcA, arcB) {
      var ia = arcA[0], ib = arcB[0],
          ja = arcA[1], jb = arcB[1],
          n = ja - ia;
      if (n !== jb - ib) return false;
      var ka = findMinimumOffset(arcA),
          kb = findMinimumOffset(arcB);
      for (var i = 0; i < n; ++i) {
        if (!equalPoint(coordinates[ia + (i + ka) % n], coordinates[ib + (i + kb) % n])) return false;
      }
      return true;
    }

    function reverseEqualRing(arcA, arcB) {
      var ia = arcA[0], ib = arcB[0],
          ja = arcA[1], jb = arcB[1],
          n = ja - ia;
      if (n !== jb - ib) return false;
      var ka = findMinimumOffset(arcA),
          kb = n - findMinimumOffset(arcB);
      for (var i = 0; i < n; ++i) {
        if (!equalPoint(coordinates[ia + (i + ka) % n], coordinates[jb - (i + kb) % n])) return false;
      }
      return true;
    }

    function findMinimumOffset(arc) {
      var start = arc[0],
          end = arc[1],
          mid = start,
          minimum = mid,
          minimumPoint = coordinates[mid];
      while (++mid < end) {
        var point = coordinates[mid];
        if (point[0] < minimumPoint[0] || point[0] === minimumPoint[0] && point[1] < minimumPoint[1]) {
          minimum = mid;
          minimumPoint = point;
        }
      }
      return minimum - start;
    }

    return topology;
  }

  function cut(topology) {
    var junctionByPoint = join(topology),
        coordinates = topology.coordinates,
        lines = topology.lines,
        rings = topology.rings;

    for (var i = 0, n = lines.length; i < n; ++i) {
      var line = lines[i],
          lineMid = line[0],
          lineEnd = line[1];
      while (++lineMid < lineEnd) {
        if (junctionByPoint.get(coordinates[lineMid])) {
          var next = {0: lineMid, 1: line[1]};
          line[1] = lineMid;
          callback({type: "cut:line", coordinates: [coordinates.slice(line[0], line[1] + 1), coordinates.slice(next[0], next[1] + 1)]});
          line = line.next = next;
        }
      }
    }

    for (var i = 0, n = rings.length; i < n; ++i) {
      var ring = rings[i],
          ringStart = ring[0],
          ringMid = ringStart,
          ringEnd = ring[1],
          ringFixed = junctionByPoint.get(coordinates[ringStart]);
      while (++ringMid < ringEnd) {
        if (junctionByPoint.get(coordinates[ringMid])) {
          if (ringFixed) {
            var next = {0: ringMid, 1: ring[1]};
            ring[1] = ringMid;
            callback({type: "cut:ring", coordinates: [coordinates.slice(ring[0], ring[1] + 1), coordinates.slice(next[0], next[1] + 1)]});
            ring = ring.next = next;
          } else {
            rotateArray(coordinates, ringStart, ringEnd, ringEnd - ringMid);
            coordinates[ringEnd] = coordinates[ringStart];
            ringFixed = true;
            callback({type: "cut:rotate", coordinates: coordinates.slice(ring[0], ring[1] + 1)});
            ringMid = ringStart;
          }
        }
      }
    }

    return topology;
  };

  function rotateArray(array, start, end, offset) {
    reverse(array, start, end);
    reverse(array, start, start + offset);
    reverse(array, start + offset, end);
  }

  function reverse(array, start, end) {
    for (var mid = start + ((end-- - start) >> 1), t; start < mid; ++start, --end) {
      t = array[start], array[start] = array[end], array[end] = t;
    }
  }

  function geomify(objects) {

    function geomifyObject(object) {
      return (object && geomifyObjectType.hasOwnProperty(object.type)
          ? geomifyObjectType[object.type]
          : geomifyGeometry)(object);
    }

    function geomifyFeature(feature) {
      var geometry = feature.geometry;
      if (geometry == null) {
        feature.type = null;
      } else {
        geomifyGeometry(geometry);
        feature.type = geometry.type;
        if (geometry.geometries) feature.geometries = geometry.geometries;
        else if (geometry.coordinates) feature.coordinates = geometry.coordinates;
      }
      delete feature.geometry;
      return feature;
    }

    function geomifyGeometry(geometry) {
      if (!geometry) return {type: null};
      if (geomifyGeometryType.hasOwnProperty(geometry.type)) geomifyGeometryType[geometry.type](geometry);
      return geometry;
    }

    var geomifyObjectType = {
      Feature: geomifyFeature,
      FeatureCollection: function(collection) {
        collection.type = "GeometryCollection";
        collection.geometries = collection.features;
        collection.features.forEach(geomifyFeature);
        delete collection.features;
        return collection;
      }
    };

    var geomifyGeometryType = {
      GeometryCollection: function(o) {
        var geometries = o.geometries, i = -1, n = geometries.length;
        while (++i < n) geometries[i] = geomifyGeometry(geometries[i]);
      },
      MultiPoint: function(o) {
        if (!o.coordinates.length) {
          o.type = null;
          delete o.coordinates;
        } else if (o.coordinates.length < 2) {
          o.type = "Point";
          o.coordinates = o.coordinates[0];
        }
      },
      LineString: function(o) {
        if (!o.coordinates.length) {
          o.type = null;
          delete o.coordinates;
        }
      },
      MultiLineString: function(o) {
        for (var lines = o.coordinates, i = 0, N = 0, n = lines.length; i < n; ++i) {
          var line = lines[i];
          if (line.length) lines[N++] = line;
        }
        if (!N) {
          o.type = null;
          delete o.coordinates;
        } else if (N < 2) {
          o.type = "LineString";
          o.coordinates = lines[0];
        } else {
          o.coordinates.length = N;
        }
      },
      Polygon: function(o) {
        for (var rings = o.coordinates, i = 0, N = 0, n = rings.length; i < n; ++i) {
          var ring = rings[i];
          if (ring.length) rings[N++] = ring;
        }
        if (!N) {
          o.type = null;
          delete o.coordinates;
        } else {
          o.coordinates.length = N;
        }
      },
      MultiPolygon: function(o) {
        for (var polygons = o.coordinates, j = 0, M = 0, m = polygons.length; j < m; ++j) {
          for (var rings = polygons[j], i = 0, N = 0, n = rings.length; i < n; ++i) {
            var ring = rings[i];
            if (ring.length) rings[N++] = ring;
          }
          if (N) {
            rings.length = N;
            polygons[M++] = rings;
          }
        }
        if (!M) {
          o.type = null;
          delete o.coordinates;
        } else if (M < 2) {
          o.type = "Polygon";
          o.coordinates = polygons[0];
        } else {
          polygons.length = M;
        }
      }
    };

    for (var key in objects) {
      objects[key] = geomifyObject(objects[key]);
    }

    return objects;
  }

  function index(objects) {
    var topology = dedup(cut(extract(objects))),
        coordinates = topology.coordinates,
        indexByArc = hashtable(topology.arcs.length, hashArc, equalArc);

    objects = topology.objects; // for garbage collection

    topology.arcs = topology.arcs.map(function(arc, i) {
      indexByArc.set(arc, i);
      return coordinates.slice(arc[0], arc[1] + 1);
    });

    delete topology.coordinates;
    coordinates = null;

    function indexGeometry(geometry) {
      if (geometry && indexGeometryType.hasOwnProperty(geometry.type)) indexGeometryType[geometry.type](geometry);
    }

    var indexGeometryType = {
      GeometryCollection: function(o) { o.geometries.forEach(indexGeometry); },
      LineString: function(o) { o.arcs = indexArcs(o.arcs); },
      MultiLineString: function(o) { o.arcs = o.arcs.map(indexArcs); },
      Polygon: function(o) { o.arcs = o.arcs.map(indexArcs); },
      MultiPolygon: function(o) { o.arcs = o.arcs.map(indexMultiArcs); }
    };

    function indexArcs(arc) {
      var indexes = [];
      do {
        var index = indexByArc.get(arc);
        indexes.push(arc[0] < arc[1] ? index : ~index);
      } while (arc = arc.next);
      return indexes;
    }

    function indexMultiArcs(arcs) {
      return arcs.map(indexArcs);
    }

    for (var key in objects) {
      indexGeometry(objects[key]);
    }

    return topology;
  }

  function hashArc(arc) {
    var i = arc[0], j = arc[1], t;
    if (j < i) t = i, i = j, j = t;
    return i + 31 * j;
  }

  function equalArc(arcA, arcB) {
    var ia = arcA[0], ja = arcA[1],
        ib = arcB[0], jb = arcB[1], t;
    if (ja < ia) t = ia, ia = ja, ja = t;
    if (jb < ib) t = ib, ib = jb, jb = t;
    return ia === ib && ja === jb;
  }

  function noop() {}

  function delta(topology) {
    var arcs = topology.arcs,
        i = -1,
        n = arcs.length;

    while (++i < n) {
      var arc = arcs[i],
          j = 0,
          m = arc.length,
          point = arc[0],
          x0 = point[0],
          y0 = point[1],
          x1,
          y1;
      while (++j < m) {
        point = arc[j];
        x1 = point[0];
        y1 = point[1];
        arc[j] = [x1 - x0, y1 - y0];
        x0 = x1;
        y0 = y1;
      }
    }

    return topology;
  }

  geomify(objects);
  var bbox = bounds(objects);
  if (Q) var transform = quantize(objects, bbox, Q);
  var topology = index(objects);
  topology.bbox = bbox;
  if (Q) topology.transform = transform, delta(topology);
  return topology;
}
