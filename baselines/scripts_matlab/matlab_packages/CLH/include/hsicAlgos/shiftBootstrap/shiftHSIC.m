function [test,defaultValues] = shiftHSIC(X,Y)
defaultValues.alpha = 0.05;

[head,tail] =  esitmateHeadAndTail(X,Y); 

defaultValues.head = head;
defaultValues.tail = tail;

defaultValues.sigx=median_heur(X);
defaultValues.sigy=median_heur(Y);

test = customShiftHSIC(X,Y, defaultValues.alpha, head,tail,...
    defaultValues.sigx,defaultValues.sigy);

end
