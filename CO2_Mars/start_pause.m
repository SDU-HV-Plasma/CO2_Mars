function start_pause
h=findobj(1,'tag','pausebutton');

if(strcmp(get(h,'string'),'Pause')) %simulation is running ->stop it
    set(h,'string','Continue');
else
    set(h,'string','Pause');
end
drawnow