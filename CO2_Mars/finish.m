function finish
h=findobj(1,'tag','endbutton');
set(h,'userdata',1);
h=findobj(1,'tag','pausebutton');
set(h,'string','Pause');