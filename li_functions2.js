function checkAllg1() {
     var checkboxes = document.getElementsByTagName('input');
	 var val = null;
     for (var i = 1; i < checkboxes.length; i=i+2) {
         if (checkboxes[i].type == 'checkbox') {
             if (val === null) val = checkboxes[i].checked;
             checkboxes[i].checked = val;
           //checkboxes[i].value = "true";  
         }
     }
    //document.getElementById('s1_12').value = "true";
 }

function checkAllg2() {
     var checkboxes = document.getElementsByTagName('input');
         var val = null;
     for (var i = 2; i < checkboxes.length; i=i+2) {
         if (checkboxes[i].type == 'checkbox') {
             if (val === null) val = checkboxes[i].checked;
             checkboxes[i].checked = val;
         }
     }
 }


