% Function to call findjobj by Yair Altman otherwise produce error
function [jObject,hObject] = javahandle()
    try
        jObject = findjobj_fast(gco);
        try %is it a multi-line edit box?
            jObject = handle(jObject.getViewport.getView, 'CallbackProperties');
        catch
        end
        hObject = gco;
    catch ME
        jObject = [];
        hObject = [];
        errordlg('Copy Paste requires findjobj by Yair Altman. Download it from Matlab File Exchange', 'Error','modal');
    end
end