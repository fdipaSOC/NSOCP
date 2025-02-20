% Function for compute: AUC, Accuracy, sensitivity, specificity, precision and F value 
% Input: 
%        Prediction - predictive  label
%        Yt         - real label
% Output:
%        AUC:   Area under curve (balance accuracy)
%        Accu:  Accuracy
%        Sens:  Sensitivity or Recall
%        Spec:  Specifity
%        cm:  Confusion Matrix
%        Prec: Precision 
%        F : F value 

function [AUC,Accu,M]=medi_auc_accu(Predict,Yt)

    Accu=1-(sum(Predict(:)~=Yt)/length(Predict(:)));
    
    tPos=sum(Yt==1 & Predict(:)==1);
    tNeg=sum(Yt==-1 & Predict(:)==-1);
    fPos=sum(Yt==-1 & Predict(:)==1);
    fNeg=sum(Yt==1 & Predict(:)==-1);
    
    Sens=tPos/(tPos+fNeg);
    Spec=tNeg/(fPos+tNeg);
    AUC=(Sens+Spec)/2;
    
    M.cm=[tNeg, fPos; fNeg, tPos];
    Prec=tPos/(tPos+fPos);
    M.F= 2*Prec*Sens/(Prec+Sens);
    M.Prec=Prec;
