const protoOXplot = { 
  width: 0.75 * width, 
  height: height, 
  margin: oxMargin,
  doses: data.doses,
  dunit: data.dunit,
  mtds: data.mtd,
  stop_esc: data.stop_esc,
  Npts: data.mtd.length,
  Nperiods: data.mtd.length/3 + 2,
};

function oxFactory(opts, proto = protoOXplot) {
  const oxPlot = Object.assign({}, proto, opts);

  const margin = oxPlot.margin;
  const horizMargins = margin.left + margin.right;
  const vertMargins = margin.top + margin.bottom;
  oxPlot.svg = div.select('#ox-plot')
      .append('svg')
      .attr('id', oxPlot.id || 'oxPlot')
      .attr('viewBox',`0 0 ${oxPlot.width+horizMargins} ${oxPlot.height+vertMargins}`)
      .attr('width', oxPlot.width + horizMargins) 
      .attr('height', oxPlot.height + vertMargins)
    .append('g')
      .attr('transform', `translate(${margin.left}, ${margin.top})`);

  return oxPlot;
}

function renderOXplot(opts) {
  const oxPlot = oxFactory(opts);

  const x = d3.scaleLinear()
    .domain([0.5, 3*oxPlot.Nperiods+0.5])
    .range([0, oxPlot.width]);

  const dsteps = oxPlot.doses.length - 1;
  
  const y_range = d3.range(oxPlot.doses.length)
    .map(d => (((dsteps+0.1)-d)*oxPlot.height/(dsteps+0.25)));
  
  const y = d3.scaleOrdinal()
    .range(y_range)
    .domain(d3.range(1, oxPlot.doses.length+1));

  const y2 = d3.scaleLog()
    .range(y_range)
    .domain(oxPlot.doses);
  
  const xAxis = d3.axisBottom().scale(x);
  xAxis.tickValues(d3.range(1,oxPlot.Npts+1));
  const yAxis = d3.axisLeft().scale(y);
  yAxis.tickValues(d3.range(1, oxPlot.doses.length+1));

  const y2Axis = d3.axisRight().scale(y2);
  y2Axis.tickValues(oxPlot.doses);
  y2Axis.tickFormat(d3.format('.2f'));

  oxPlot.svg.append('g')
    .attr('class','axis participant')
    .attr('transform', `translate(0, ${oxPlot.height})`)
    .call(xAxis);

  oxPlot.svg.append('g')
    .attr('class','axis doselevel')
    .attr('transform', 'translate(0, 0)')
    .call(yAxis);

  oxPlot.svg.append('g')
    .attr('class','axis doseabsolute')
    .attr('transform', `translate(${oxPlot.width}, 0)`)
    .call(y2Axis);

  var tooltip = div.append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

  // Thanks for the following tooltip code go to
  // http://www.d3noob.org/2013/01/adding-tooltips-to-d3js-graph.html
//  d3.csv('data/mtd.csv', data => {
    // TODO: Understand the scoping rules that have required me
    //       to attach mtds and dunit to oxPlot object, instead
    //       of having these available through -data-.
    const mtd = function(id){
      const mtd_i = oxPlot.mtds.filter(d => d.id == id)[0].mtd;
      return d3.format('.3f')(mtd_i);
    };
    oxPlot.svg.selectAll('.axis.participant .tick text')
        .on('mouseover', d => {
            showSeries(d);
            tooltip.transition()
                .duration(200)
                .style("opacity", 0.9);
            tooltip.html(`<b>MTD<sub>(i=${d})</sub> =` + mtd(d) + '</b>')
                .style("left", (d3.event.pageX - oxPlot.margin.left + 5) + "px")
                .style("top", (d3.event.pageY - oxPlot.margin.top + 20) + "px");
          })
        .on('mouseout', d => {
            unshowSeries(d);
            tooltip.transition()
                .duration(500)
                .style("opacity", 0);
          })
      .append('div')
        .attr('title','MTDi = ---');
//  });

  const topAxis = d3.axisBottom().scale(
    d3.scaleLinear()
      .domain([0.5, oxPlot.Nperiods+0.5])
      .range([0, oxPlot.width])
  );
  topAxis.tickValues(d3.range(1,oxPlot.Nperiods+1));

  oxPlot.svg.append('g')
    .attr('class','axis')
    .attr('transform', 'translate(0, -20)')
    .call(topAxis);

  const mouseoverPeriodEnd = function(per) {
     oxPlot.svg.selectAll('line.period-end[period="'+per+'"]')
       .classed('salient', true);
  };
  const mouseoutPeriodEnd = function(per) {
     oxPlot.svg.selectAll('line.period-end[period="'+per+'"]')
       .classed('salient', false);
  };

  // Draw the DLT assessment period boundaries
  // TODO: Consider a redo using inner axis ticks
  const periodlines = d3.range(oxPlot.Nperiods).map(i => ({x:(3.5 + 3*i), per:(i+1)}));
  oxPlot.svg.append('g').selectAll('.period-end')
      .data(periodlines)
    .enter().append('line')
      .attr('class','period-end')
      .attr('period', d => d.per)
      .attr('x1', d => x(d.x))
      .attr('x2', d => x(d.x))
      .attr('y1', y(1))
      .attr('y2', y(oxPlot.doses.length))
      // Use mouseover/out to intimate that period lines are clickable
      .on('mouseover', d => { mouseoverPeriodEnd(d.per); })
      .on('mouseout', d => { mouseoutPeriodEnd(d.per); })
      .on('click', d => {
          // Dim OX plot to right of clicked period-end marker
          oxPlot.svg.selectAll('.dosemarker')
            .filter(e => +e.per <= d.per)
            .style('stroke-opacity', 0.7);
          oxPlot.svg.selectAll('.dosemarker')
            .filter(e => +e.per > d.per)
            .style('stroke-opacity', 0.2);
          // Dim any 'future' stop_esc line to avoid
          // false suggestion it was pre-planned: 
          oxPlot.svg.selectAll('.stop-esc')
            .filter(e => +e.per <= d.per)
            .style('stroke-opacity', 1.0);
          oxPlot.svg.selectAll('.stop-esc')
            .filter(e => +e.per > d.per)
            .style('stroke-opacity', 0.3);
          // Select (upper, surv, lower) for the clicked period
          dsPlot.svg.selectAll('.ds-line path')
            .style('visibility','hidden');
          dsPlot.svg.selectAll('.ds-line path[period="'+d.per+'"]')
            .style('visibility','visible');
      });

  // Draw the stop_esc indicator(s), of which there may be 0, 1 or several.
  // Initially, I employ the same dashed-line period divider overlay as in
  // the original OXplot based on lattice graphics.
  // TODO: Consider a (more intuitive?) horizontal escalation 'barrier'.
  if(Object.keys(oxPlot.stop_esc).length === 0) oxPlot.stop_esc = [];
  const esc_stops = (oxPlot.stop_esc).map(i => ({x:(0.5 + 3*i), per:i}));
  oxPlot.svg.append('g').selectAll('.stop-esc')
      .data(esc_stops)
    .enter().append('line')
      .attr('class','stop-esc')
      .attr('period', d => d.per)
      .attr('x1', d => x(d.x))
      .attr('x2', d => x(d.x))
      .attr('y1', y(1))
      .attr('y2', y(oxPlot.doses.length));

  const showSeries = function(pid) {
    oxPlot.svg.select('g.axis').selectAll('g.tick')
        .filter(d => d==`${pid}`)
        .style('font-weight', 'bold')
        .style('font-size', 14);
    oxPlot.svg.selectAll('.dosemarker[participant="'+pid+'"]')
      .classed('salient', true);
    oxPlot.svg.selectAll('.trace[participant="'+pid+'"]')
      .classed('salient', true);
    //dsPlot.svg.selectAll('.ds-line path')
    //  .style('visibility','hidden');
    dsPlot.svg.selectAll('.ds-pointer path[participant="'+pid+'"]')
      .style('visibility', 'visible');
  };

  const unshowSeries = function(pid) {
    oxPlot.svg.select('g.axis').selectAll('g.tick')
        .style('font-weight', 'normal')
      .style('font-size', 10);
    oxPlot.svg.selectAll('.dosemarker[participant="'+pid+'"]')
      .classed('salient', false);
    oxPlot.svg.selectAll('.trace[participant="'+pid+'"]')
      .classed('salient', false);
    dsPlot.svg.selectAll('.ds-pointer path[participant="'+pid+'"]')
      .style('visibility', 'hidden');
  };

  const average = arr => arr.reduce( ( p, c ) => p + c, 0 ) / arr.length;

  const ptMarker = d3.symbol()
    .type(d => (d.dlt ? symbolExx : d3.symbolCircle))
    .size(height/2); // TODO: Isn't this -size- parameter awfully large?

  const trace = d3.line()
      .x(d => x((d.id-1)%3 + 3*(d.period-1) + 1))
      .y(d => y(d.dose));

  const takeUnique = function(v, i, a) { return a.indexOf(v) === i };

//  d3.csv('data/trial.csv', data => {
    // -jitters- maps (period, dose) keys to arrays of cohort numbers
    const jitters = d3.nest().key(d => d.period).key(d => d.dose)
      .rollup(leaves => leaves.map(v => cohort(v.id)).filter(takeUnique))
      .entries(data.trial);
    //console.log(jitters);
    //console.log(jitters.filter(d => d.key==3)[0].values
    //                   .filter(d => d.key==2)[0].value);
    // Draw the participant-wise traces
    const pidwise = d3.nest().key(d => d.id).entries(data.trial);
    pidwise.forEach(i => {
      oxPlot.svg.append('g')
        .append('path')
        .attr('participant', i.key)
        .attr('class', 'trace')
        .attr('d', trace(i.values));
    });
    // Draw the dosemarkers (i.e., O's and X's)
    data.trial.forEach(i => {
      // -pdcohs- is to be an array of cohort numbers represented at (i.period, i.dose)
      const pdcohs = jitters.filter(d => d.key==i.period)[0].values
                            .filter(d => d.key==i.dose)[0].value;
      var spread = d3.range(0, pdcohs.length);
      spread = spread.map(d => d - d3.mean(spread));
      const yOffset = 10*spread[pdcohs.indexOf(cohort(i.id))];
      var abscissa = x((i.id-1)%3 + 3*(i.period-1) + 1);
      var ordinate = y(+i.dose) + yOffset;
      oxPlot.svg.append('path')
          .datum({per: `${i.period}`}) // TODO: Use .data() idiom
          .attr('class', 'dosemarker')
          .attr('d', ptMarker(i))
          .attr('stroke', colorForID(i.id))
          .attr('transform', `translate(${abscissa}, ${ordinate})`)
          .attr('period', `${i.period}`)
          .attr('participant', `${i.id}`)
          .attr('dose-level', `${i.dose}`)
          .attr('dlt-seen', `${i.dlt}`)
          .on('mouseover', () => { showSeries(i.id); })
          .on('mouseout', () => { unshowSeries(i.id); });
    });
//  }); 

  oxPlot.svg.append("text")
    .attr("class", "x label")
    .attr("text-anchor", "middle")
    .attr("x", x(15))
    .attr("y", height + 35)
    .text("Participant number");

  oxPlot.svg.append("text")
    .attr("class", "y label")
    .attr("text-anchor", "middle")
    .attr("x", -height*0.5)
    .attr("y", -35)
    .attr("dy", ".75em")
    .attr("transform", "rotate(-90)")
    .text("Dose level");

  oxPlot.svg.append("text")
    .attr("class", "per label")
    .attr("text-anchor", "middle")
    .attr("x", x(15))
    .attr("y", 10)
    .text("DLT assessment period");

  oxPlot.svg.append("text")
    .attr("class", "y label")
    .attr("text-anchor", "middle")
    .attr("x", -height*0.5)
    .attr("y", oxPlot.width + 35)
    .attr("dy", ".75em")
    .attr("transform", "rotate(-90)")
    .text("Dose (" + data.dunit + ")");

  return oxPlot;
}
