const oxFract = 0.75;
const dsFract = 1.0 - oxFract;

const oxds_divs = [
  {id:"ox-plot", class:"svg-container"},
  {id:"ds-plot", class:"svg-container"}];

// Provide DIVs for the OX and DS plots:
div.selectAll("div").data(oxds_divs).enter().append("div")
    .attr("id", d => d.id)
    .attr("class", d => d.class);

// TODO: Should I use 'classed()' above, instead?

div.selectAll("div").selectAll("svg").remove(); // clean slate enables update!
    
div.select('#ox-plot')
    .style('width', oxFract * width + horiz(oxMargin) + 'px')
    .style('height', height + vert(oxMargin) + 'px')
    .style('position','absolute')
    .style('top', 0+'px')
    .style('left', 0+'px');

div.select('#ds-plot')
    .style('width', dsFract * width + horiz(dsMargin) + 'px')
    .style('height', height + vert(dsMargin) + 'px')
    .style('position','absolute')
    .style('top', 0+'px')
    .style('left', oxFract * width + horiz(oxMargin) + 'px');

// Establish some shared infrastructure to unify OX & DS
// TODO: Factor lots more 'infrastructure' up here, e.g. shared axes.
const cohort = id => Math.floor((+id - 1) / 3) + 1;
// Colors from brewer.pal(4,"Dark2")
const colorCycle = ['#1B9E77','#D95F02','#7570B3','#E7298A'];
const colorForID = id => colorCycle[(cohort(id)-1) % colorCycle.length];

// TODO: Consider passing the appropriate DIVs to these 'render*' functions.
//       This would insulate these functions from need for esoteric knowledge
//       of the -div- variable!
// OTOH: Maybe it's better not to pretend I'm achieving more modularity than
//       this programming environment truly supports. That is, perhaps it's
//       best to acknowledge the prevailing idioms, and treat r2d3-provided
//       variables as globals. I can expect to gain little from a superficial
//       simulacrum of modularity!
this.oxPlot = renderOXplot({width: oxFract*width,
                              margin: oxMargin,
                              });
this.dsPlot = renderDSplot({width: dsFract*width,
                              margin: dsMargin,
                              });

div.select('.footer')
    .style('top', height + vert(oxMargin) +'px');