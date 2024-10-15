function [  ] = celledit( hobj,evnt,uf,paraanz  )
% Special cell edit to enable blockwise (de)activation of active paramters

modifiers = get(uf,'CurrentModifier');
CtrlPressed  = ismember('control', modifiers);
act=get(hobj,'data');
actProp=act.Properties.VariableNames(evnt.Indices(2));

if CtrlPressed && strcmp(actProp{1},'Active')
    
    act(find(paraanz==paraanz(evnt.Indices(1))),evnt.Indices(2))={evnt.NewData}; %#ok<FNDSB>
else
    act(evnt.Indices(1),evnt.Indices(2))={evnt.NewData};
end
set(hobj,'data',act)
end

