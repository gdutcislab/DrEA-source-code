function mop = testmop( testname )
    path('test',path);
    mop=get_structure( 'testmop' );
    switch testname
        case 'WFG1'
            mop=wfg1(mop);
        case 'WFG2'
            mop=wfg2(mop);
        case 'WFG3'
            mop=wfg3(mop);
        case 'WFG4'
            mop=wfg4(mop);
        case 'WFG5'
            mop=wfg5(mop);
        case 'WFG6'
            mop=wfg6(mop);
        case 'WFG7'
            mop=wfg7(mop);
        case 'WFG8'
            mop=wfg8(mop);
        case 'WFG9'
            mop=wfg9(mop);
        case 'DTLZ1'
            mop=DTLZ1(mop);
        case 'DTLZ2'
            mop=DTLZ2(mop);
        case 'DTLZ3'
            mop=DTLZ3(mop);
        case 'DTLZ4'
            mop=UF4(mop);
        case 'UF5'
            mop=UF5(mop);
        case 'UF6'
            mop=UF1(mop);
        case 'UF6'
            mop=UF7(mop);
        case 'UF7'
            mop=UF8(mop);
        case 'UF8'
            mop=UF9(mop);
        case 'UF9'
            mop=UF1(mop);
        case 'UF10'
            mop=UF10(mop);            
        case 'MOP2'
            mop=MOP2(mop);
        case 'MOP3'
            mop=MOP3(mop);
        case 'MOP4'
            mop=MOP4(mop);
        case 'MOP5'
            mop=MOP5(mop);
        case 'MOP6'
            mop=MOP6(mop);
        case 'MOP7'
            mop=MOP7(mop);
        case 'ZDT1'
            mop=ZDT1(mop);
        case 'ZDT2'
            mop=ZDT2(mop);
        otherwise
            error('Undefined test problem name');
    end     
end
