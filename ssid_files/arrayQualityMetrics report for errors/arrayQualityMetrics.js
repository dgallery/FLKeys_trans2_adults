// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, true, false, false, false, false, false, false, false, true, false, false, false, false, false, false, true, false, false, true, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "SA10D", "SA10D", "deep", "S1" ], [ "2", "SA10N", "SA10N", "nearshore", "S1" ], [ "3", "SA10O", "SA10O", "offshore", "S1" ], [ "4", "SA11D", "SA11D", "deep", "S1" ], [ "5", "SA11O", "SA11O", "offshore", "S1" ], [ "6", "SA12D", "SA12D", "deep", "S1" ], [ "7", "SA12N", "SA12N", "nearshore", "S1" ], [ "8", "SA12O", "SA12O", "offshore", "S1" ], [ "9", "SA13D", "SA13D", "deep", "S1" ], [ "10", "SA13N", "SA13N", "nearshore", "S1" ], [ "11", "SA14D", "SA14D", "deep", "D1" ], [ "12", "SA14N", "SA14N", "nearshore", "S2" ], [ "13", "SA14O", "SA14O", "offshore", "S1" ], [ "14", "SA15D", "SA15D", "deep", "D1" ], [ "15", "SA15N", "SA15N", "nearshore", "S1" ], [ "16", "SA15O", "SA15O", "offshore", "S1" ], [ "17", "SA16D", "SA16D", "deep", "S1" ], [ "18", "SA16N", "SA16N", "nearshore", "S1" ], [ "19", "SA16O", "SA16O", "offshore", "S1" ], [ "20", "SA17D", "SA17D", "deep", "D1" ], [ "21", "SA17N", "SA17N", "nearshore", "S1" ], [ "22", "SA17O", "SA17O", "offshore", "S1" ], [ "23", "SA18D", "SA18D", "deep", "S1" ], [ "24", "SA18N", "SA18N", "nearshore", "S1" ], [ "25", "SA18O", "SA18O", "offshore", "S1" ], [ "26", "SA19D", "SA19D", "deep", "S1" ], [ "27", "SA19N", "SA19N", "nearshore", "S1" ], [ "28", "SA19O", "SA19O", "offshore", "S1" ], [ "29", "SA1D", "SA1D", "deep", "S1" ], [ "30", "SA1N", "SA1N", "nearshore", "S1" ], [ "31", "SA20D", "SA20D", "deep", "S1" ], [ "32", "SA20N", "SA20N", "nearshore", "S1" ], [ "33", "SA20O", "SA20O", "offshore", "S1" ], [ "34", "SA2D", "SA2D", "deep", "D1" ], [ "35", "SA2N", "SA2N", "nearshore", "S1" ], [ "36", "SA2O", "SA2O", "offshore", "S1" ], [ "37", "SA3D", "SA3D", "deep", "D2" ], [ "38", "SA3N", "SA3N", "nearshore", "S1" ], [ "39", "SA3O", "SA3O", "offshore", "S1" ], [ "40", "SA4D", "SA4D", "deep", "D1" ], [ "41", "SA4N", "SA4N", "nearshore", "S1" ], [ "42", "SA4O", "SA4O", "offshore", "S1" ], [ "43", "SA5D", "SA5D", "deep", "D1" ], [ "44", "SA5N", "SA5N", "nearshore", "S1" ], [ "45", "SA5O", "SA5O", "offshore", "S1" ], [ "46", "SA6D", "SA6D", "deep", "S1" ], [ "47", "SA6N", "SA6N", "nearshore", "S1" ], [ "48", "SA6O", "SA6O", "offshore", "S1" ], [ "49", "SA7D", "SA7D", "deep", "D2" ], [ "50", "SA7N", "SA7N", "nearshore", "S1" ], [ "51", "SA7O", "SA7O", "offshore", "S1" ], [ "52", "SA8D", "SA8D", "deep", "S1" ], [ "53", "SA8N", "SA8N", "nearshore", "S1" ], [ "54", "SA8O", "SA8O", "offshore", "S1" ], [ "55", "SA9D", "SA9D", "deep", "S1" ], [ "56", "SA9N", "SA9N", "nearshore", "S1" ], [ "57", "SA9O", "SA9O", "offshore", "S1" ] ];
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
