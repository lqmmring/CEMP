function [RS,Jc,FM,F1]=evalution(X,idx, label)
    Row = size(X,1);
    ss=0;sd=0;ds=0;dd=0;
    for i=1:Row-1
        for j=i+1:Row
            if idx(i)==idx(j) && label(i)==label(j)
                ss=ss+1;
            elseif idx(i)==idx(j) && label(i)~=label(j)
                sd=sd+1;
            elseif idx(i)~=idx(j) && label(i)==label(j)
                ds=ds+1;
            else
                dd=dd+1;
            end
        end
    end
    trn=ss+sd+ds+dd;
    pr=ss/(ss+sd)+ss/(ss+ds);
%     fprintf('总关系数目：%d   \n',trn);
    RS=(ss+dd)/trn;
    Jc=ss/(ss+sd+ds);
    FM=sqrt((ss^2)/((ss+sd)*(ss+ds)));
    F1=2*(ss/(ss+sd))*(ss/(ss+ds))/pr;
%     fprintf('Rand指标：%f   \n',RS);
%     fprintf('Jaccard指标：%f   \n',Jc);
%     fprintf('FM指标：%f   \n',FM);
%     fprintf('F1指标：%f   \n',F1);
end