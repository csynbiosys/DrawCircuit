function model = gen_promoter_cds_problem(promoters,CDSs,RBSs,expData)

N_PROMOTERS=length(promoters);
N_CDS=length(CDSs);

eqns={};
promEqns={};
parNames={};
parValues=[];
states={};

for i=1:length(promoters)
    
    eval(promoters{i});
    promEqns={promEqns{:} [promoters{i} '=' transcriptionFunction ';']};
    parNames={parNames{:},parsNamesTranscription{:}};
    parValues=[parValues parsValsTranslation];
    
end

promEqns=promEqns';

counter=0;
for i=1:length(CDSs)
    
    eval(CDSs{i});
    eval(RBSs{i})
    eqn=['d' CDSs{i} '='  ] ;
    
    states={states{:} CDSs{i}};
    
    parNames={parNames{:},['translationEffiency_' promoters{i}]};
    parValues=[parValues translationEffiency];
    
    for j=1:length(promoters)
        eqn=sprintf('%s\n',eqn);
        counter=counter+1;
        eqn =[eqn sprintf('\t%s',[ '+ translationEffiency_' promoters{i} '* Y' num2str(counter) ' *' promoters{i}])];
    end
    
    eqn=[eqn ';'];
    eqns={eqns{:} eqn};
    
end

counter=0;
for i=1:length(CDSs)
    
    eval(CDSs{i});
    
    for j=1:length(codes)
        
        if strcmp(codes{j},'protein')
            
            eqn=['d' names{j} '='  translationParameter '*' CDSs{i}] ;
            states={names{j} CDSs{i}};
            eval(names{j});
            parNames={parNames{:} parNamesCDS{:}};
            parValues=[parValues parValuesCDS];
            eqn=[eqn '-' parName '*' transcriptDegradationParameter];
               
        else
            
            stop('Invalid sequece type!');
            
        end
        
    end
    
    eqn=[eqn ';'];
    eqns={eqns{:} eqn};
end

NBINS=counter;
NCONT=length(parValues);
parValues=[parValues zeros(1,NBINS)];

for i=1:NBINS
    parNames={parNames{:} ['Y' numstr(i)]};
end


eqns=eqns';
eqns={promEqns{:} eqns{:}};

model.parNames=parNames;
model.parValues=parValues;
model.eqns=eqns;
model.states=states;

end
