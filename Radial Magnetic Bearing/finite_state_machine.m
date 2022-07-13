y=TRIP(0)
function y = TRIP(u)
S1=0;
S2=1;
persistent current_state;
if isempty(current_state)
    current_state=S1;
end

switch(current_state)
    case S1,
        if(u)
            y=false;
            current_state=S1;
        else
            y=true;
            current_state=S2;
        end
        
    case S2,
        if(u)
            y=false;
            current_state=S2;
        else
            y=true;
            current_state=S2;
        end
    otherwise,
        y=false;
end
if current_state==S1
    y=1;
else
    current_state==S2
    y=0;
end

if(current_state==1)
    y=25;
else
    y=12
end
end