// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, true, false, false, false ];
var arrayMetadata    = [ [ "1", "MA10D", "MA10D", "deep", "D3" ], [ "2", "MA10N", "MA10N", "nearshore", "O2" ], [ "3", "MA10O", "MA10O", "offshore", "N1" ], [ "4", "MA11N", "MA11N", "nearshore", "N1" ], [ "5", "MA12N", "MA12N", "nearshore", "N1" ], [ "6", "MA13N", "MA13N", "nearshore", "N1" ], [ "7", "MA13O", "MA13O", "offshore", "N1" ], [ "8", "MA14N", "MA14N", "nearshore", "N1" ], [ "9", "MA15D", "MA15D", "deep", "D3" ], [ "10", "MA15N", "MA15N", "nearshore", "N1" ], [ "11", "MA15O", "MA15O", "offshore", "O2" ], [ "12", "MA16N", "MA16N", "nearshore", "O2" ], [ "13", "MA16O", "MA16O", "offshore", "O2" ], [ "14", "MA17D", "MA17D", "deep", "D3" ], [ "15", "MA17O", "MA17O", "offshore", "O2" ], [ "16", "MA18N", "MA18N", "nearshore", "N1" ], [ "17", "MA18O", "MA18O", "offshore", "O2" ], [ "18", "MA19D", "MA19D", "deep", "D3" ], [ "19", "MA19N", "MA19N", "nearshore", "N1" ], [ "20", "MA1N", "MA1N", "nearshore", "O2" ], [ "21", "MA1O", "MA1O", "offshore", "N1" ], [ "22", "MA20D", "MA20D", "deep", "D3" ], [ "23", "MA2D", "MA2D", "deep", "D3" ], [ "24", "MA2N", "MA2N", "nearshore", "N1" ], [ "25", "MA3D", "MA3D", "deep", "D3" ], [ "26", "MA3O", "MA3O", "offshore", "O2" ], [ "27", "MA4D", "MA4D", "deep", "D3" ], [ "28", "MA6N", "MA6N", "nearshore", "N1" ], [ "29", "MA6O", "MA6O", "offshore", "N1" ], [ "30", "MA7O", "MA7O", "offshore", "O2" ], [ "31", "MA8D", "MA8D", "deep", "D4" ], [ "32", "MA8N", "MA8N", "nearshore", "O2" ], [ "33", "MA8O", "MA8O", "offshore", "O2" ], [ "34", "MA9N", "MA9N", "nearshore", "O2" ], [ "35", "MA9O", "MA9O", "offshore", "D3" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
