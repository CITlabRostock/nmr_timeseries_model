%% Small auxiliary function, closes GUI
function my_closereq(src,callbackdata,uf,uit,tf)
try
    setappdata(tf,'tablecontent',uit.Data)
    delete(uf)
catch
    disp('Error while close Parameter Overview')
    delete(uf)
end
end