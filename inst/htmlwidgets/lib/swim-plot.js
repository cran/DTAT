const protoSWIMplot = { 
  width: 0.75 * width, 
  height: height * 0.45, 
  margin: oxMargin,
  doses: data.doses,
  dunit: data.dunit,
  mtds: data.mtd,
  stop_esc: data.stop_esc,
  Npts: data.mtd.length,
  Nperiods: data.mtd.length/3 + 2,
  trial: data.trial,
};

function swimFactory(opts, proto = protoSWIMplot) {
  const swimPlot = Object.assign({}, proto, opts);

  const margin = swimPlot.margin;
  margin.top = 0;
  const horizMargins = margin.left + margin.right;
  const vertMargins = margin.top + margin.bottom;
  swimPlot.svg = div.select('#ox-plot')
      .append('svg')
      .attr('id', swimPlot.id || 'swimPlot')
      .attr('viewBox',`0 0 ${swimPlot.width+horizMargins} ${swimPlot.height+vertMargins}`)
      .attr('width', swimPlot.width + horizMargins) 
      .attr('height', swimPlot.height + vertMargins)
    .append('g')
      .attr('transform', `translate(${margin.left}, ${margin.top})`);

  return swimPlot;
}

/**
 * TODO:
 *v1. Add a 'weeks' y axis
 *v2. Render a single, light-blue bar for each participant & dose-period
 *v3. Make the bar intensity (alpha?) a function of period dose
 *v4. Refine sizing and margins (maybe double-subtract?)
 *v5. Extend the 'data' object beyond the 'exit rule'
 * - The key here (no pun intended!) seems to be making appropriate use
 *   of the map/filter/reduce functionality in Javascript. I get a fine
 *   opportunity at last to learn this API!
 * - First step would seem to be creating a list of {id, period, dose}
 *   by means of a *reduce* operation that obtains for me the very last
 *   recorded period for each id.
 * - For each id, I should filter the array, sort by period (descending)
 *   and take the first element of the result. From this, I should then
 *   construct an array of the censored dose periods, and append this
 *   array to the trials.
 * - Consider this code, for example [which works in console]
 *   > swimPlot.trial.filter(d => d.id == 1).sort(d => -d.period)[0]
 * 6. Add data-driven clinical event markers with characteristic symbols
 * - Probably these data can be provided in future via additional fields
 *   of the 'data' object, accessed by 'data.response' or 'data.events'.
 *   But to get started, let me just hard-code such data here in these
 *   Javascript files.
 * 7. Add a KEY in the lower, right-hand corner
*/

function renderSWIMplot(opts) {
  const swimPlot = swimFactory(opts);
  
  // TODO: Consider moving the following inside swimFactory
  const furtherDosing = function(dp){
    // dp is a 'doseperiod' of the form {id:1, period:7, dose:6, dlt:false}
    furtherDose = dp.dlt ? dp.dose-1 : dp.dose;
    if (furtherDose > 0) {
      return d3.range(dp.period, swimPlot.Nperiods).map(
        per => ({id: dp.id, period: per+1, dose: furtherDose}));
    } else {
      return [];
    }
  };
  
  var lastdose = [];
  d3.nest().key(d => d.id).entries(swimPlot.trial).forEach(i => {
    lastdose.push(furtherDosing(i.values.sort(
      (d1,d2) => d2.period - d1.period)[0]));
  });
  swimPlot.lastdose = lastdose;
  swimPlot.trial = swimPlot.trial.concat(d3.merge(lastdose));

  swimPlot.responses = [
    {id: 3, week: 8, dose: 2, response: "PR", ott: 1.3},
    {id: 6, week: 10, dose: 3, response: "CR", ott: 1.5},
    {id: 7, week: 16, dose: 4, response: "PR", ott: 2.5},
    {id: 12, week: 6, dose: 2, response: "PR", ott: 0.5},
    {id: 14, week: 16, dose: 5, response: "PR", ott: 0.5},
    {id: 18, week: 14, dose: 4, response: "CR", ott: 3.2},
    ];

  const x = d3.scaleLinear()
    .domain([0.5, 3*swimPlot.Nperiods+0.5])
    .range([0, oxPlot.width]);

  const dsteps = swimPlot.doses.length - 1;
  
  // TODO: Change y axis to *time*
  const y_range = d3.range(swimPlot.doses.length)
    .map(d => (((dsteps+0.1)-d)*swimPlot.height/(dsteps+0.25)));
  
  const y = d3.scaleLinear()
    .domain([0.0, swimPlot.Nperiods*4])
    .range([0, swimPlot.height]);

  const yAxis = d3.axisLeft().scale(y);
  yAxis.tickValues(d3.range(4, (1+swimPlot.Nperiods)*4, 4));

  swimPlot.svg.append('g')
    .attr('class','axis doselevel')
    .attr('transform', 'translate(0, 0)')
    .call(yAxis);

  var tooltip = div.append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

  const mouseoverPeriodEnd = function(per) {
     swimPlot.svg.selectAll('line.period-end[period="'+per+'"]')
       .classed('salient', true);
  };
  const mouseoutPeriodEnd = function(per) {
     swimPlot.svg.selectAll('line.period-end[period="'+per+'"]')
       .classed('salient', false);
  };

  const showSeries = function(pid) {
    swimPlot.svg.select('g.axis').selectAll('g.tick')
        .filter(d => d==`${pid}`)
        .style('font-weight', 'bold')
        .style('font-size', 14);
    swimPlot.svg.selectAll('.dosemarker[participant="'+pid+'"]')
      .classed('salient', true);
    swimPlot.svg.selectAll('.trace[participant="'+pid+'"]')
      .classed('salient', true);
    //dsPlot.svg.selectAll('.ds-line path')
    //  .style('visibility','hidden');
    dsPlot.svg.selectAll('.ds-pointer path[participant="'+pid+'"]')
      .style('visibility', 'visible');
  };

  const unshowSeries = function(pid) {
    swimPlot.svg.select('g.axis').selectAll('g.tick')
        .style('font-weight', 'normal')
      .style('font-size', 10);
    swimPlot.svg.selectAll('.dosemarker[participant="'+pid+'"]')
      .classed('salient', false);
    swimPlot.svg.selectAll('.trace[participant="'+pid+'"]')
      .classed('salient', false);
    dsPlot.svg.selectAll('.ds-pointer path[participant="'+pid+'"]')
      .style('visibility', 'hidden');
  };

  const average = arr => arr.reduce( ( p, c ) => p + c, 0 ) / arr.length;

  const trace = d3.line()
      .x(d => x((d.id-1)%3 + 3*(d.period-1) + 1))
      .y(d => y(d.dose));

  const takeUnique = function(v, i, a) { return a.indexOf(v) === i };

    // Draw the dose periods, with intensity reflecting dose
    swimPlot.trial.forEach(i => {
      swimPlot.svg.append('rect')
          .datum({per: `${i.period}`}) // TODO: Use .data() idiom
          .attr("class", "doseperiod")
          .attr("dose", i.dose)
          .attr("x", x(i.id - 0.3))
          .attr("y", 4*y(i.period - cohort(i.id)))
          .attr("width", 0.6*(x(2) - x(1)))
          .attr("height", 4*(y(2) - y(1)))
          .attr("fill", "rgba(20,50,255,"+ 0.75*i.dose/swimPlot.doses.length +")");
    });
    
    swimPlot.responses.forEach(i => {
      swimPlot.svg.append('path') // TODO: vary symbol by response
          .datum({id: i.id, week: i.week, response: i.response})
          .attr("class", "response")
          .attr('d', responseMarker().size(x(3)-x(1)))
          .attr('stroke', "black")
          .attr('transform', `translate(${x(i.id)}, ${y(i.week)})`)
          .attr("assessment", i.response)
          .attr("week", i.week)
          .attr("x", x(i.id))
          .attr("y", y(i.week - 4*cohort(i.id)));
    });

  swimPlot.svg.append("text")
    .attr("class", "x label")
    .attr("text-anchor", "middle")
    .attr("x", x(15))
    .attr("y", height + 35)
    .text("Participant number");

  swimPlot.svg.append("text")
    .attr("class", "y label")
    .attr("text-anchor", "middle")
    .attr("x", -swimPlot.height*0.5)
    .attr("y", -35)
    .attr("dy", ".75em")
    .attr("transform", "rotate(-90)")
    .text("Weeks on study");

  // Add a legend
  var legend = swimPlot.svg.append("g")
    .attr('class','events-legend')
    .attr('transform', `translate(${0.8*swimPlot.width}, ${0.8*swimPlot.height})`);
    
  legend.append('path')
    .datum({response: "CR"})
    .attr('d', responseMarker().size(x(3)-x(1)))
    .attr('stroke', "black")
    .attr('transform', "translate(5,5)");
  legend.append('text')
    .attr("x", 15)
    .attr("y", 5)
    .attr("dy", "0.3em")
    .text("CR");
  legend.append('path')
    .datum({response: "PR"})
    .attr('d', responseMarker().size(x(3)-x(1)))
    .attr('stroke', "black")
    .attr('transform', "translate(5,25)");
  legend.append('text')
    .attr("x", 15)
    .attr("y", 25)
    .attr("dy", "0.3em")
    .text("PR");


  return swimPlot;
}
