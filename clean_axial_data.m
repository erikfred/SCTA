% clean_axial_data.m
%
% Simple script for truncating and cleaning tilt and pressure data prior to
% quantifying deflation events
%

ttrim=[16000 140000 1 1 140000];
ptrim=[1 140000 1 1 140000];
if eventnum==1 % no ASHES data
    AXCC1.t_time=AXCC1.t_time(ttrim(eventnum):end);
    AXCC1.LAX=AXCC1.LAX(ttrim(eventnum):end);
    AXCC1.LAY=AXCC1.LAY(ttrim(eventnum):end);
    AXCC1.p_time=AXCC1.p_time(ptrim(eventnum):end);
    AXCC1.BDO=AXCC1.BDO(ptrim(eventnum):end);
    
    AXCC1.t_time(268645:268705)=[];AXCC1.LAX(268645:268705)=[];AXCC1.LAY(268645:268705)=[]; % jump here
    AXCC1.LAX(268645:end)=AXCC1.LAX(268645:end)-(AXCC1.LAX(268645)-AXCC1.LAX(268644));
    AXCC1.LAY(268645:end)=AXCC1.LAY(268645:end)-(AXCC1.LAY(268645)-AXCC1.LAY(268644));
    
    AXCC1.t_time(251270:251330)=[];AXCC1.LAX(251270:251330)=[];AXCC1.LAY(251270:251330)=[]; % jump here
    AXCC1.LAX(251270:end)=AXCC1.LAX(251270:end)-(AXCC1.LAX(251270)-AXCC1.LAX(251269));
    AXCC1.LAY(251270:end)=AXCC1.LAY(251270:end)-(AXCC1.LAY(251270)-AXCC1.LAY(251269));
    
    AXCC1.t_time(232230:232430)=[];AXCC1.LAX(232230:232430)=[];AXCC1.LAY(232230:232430)=[]; % jump here
    AXCC1.LAX(232230:end)=AXCC1.LAX(232230:end)-(AXCC1.LAX(232230)-AXCC1.LAX(232229));
    AXCC1.LAY(232230:end)=AXCC1.LAY(232230:end)-(AXCC1.LAY(232230)-AXCC1.LAY(232229));
    
    AXEC2.LAY=AXEC2.LAY(159000:end);
    AXEC2.t_time=AXEC2.t_time(end-length(AXEC2.LAY)+1:end);
    AXEC2.LAX=AXEC2.LAX(end-length(AXEC2.LAY)+1:end);
    AXEC2.p_time=AXEC2.p_time(ptrim(eventnum):end);
    AXEC2.BDO=AXEC2.BDO(ptrim(eventnum):end);
    
    AXEC2.t_time(59210:59270)=[];AXEC2.LAX(59210:59270)=[];AXEC2.LAY(59210:59270)=[]; % jump here
    AXEC2.LAX(59210:end)=AXEC2.LAX(59210:end)-(AXEC2.LAX(59210)-AXEC2.LAX(59209));
    AXEC2.LAY(59210:end)=AXEC2.LAY(59210:end)-(AXEC2.LAY(59210)-AXEC2.LAY(59209));
    
    AXEC2.t_time(54255:54315)=[];AXEC2.LAX(54255:54315)=[];AXEC2.LAY(54255:54315)=[]; % jump here
    AXEC2.LAX(54255:end)=AXEC2.LAX(54255:end)-(AXEC2.LAX(54255)-AXEC2.LAX(54254));
    AXEC2.LAY(54255:end)=AXEC2.LAY(54255:end)-(AXEC2.LAY(54255)-AXEC2.LAY(54254));
    
    AXID1.t_time=AXID1.t_time(ttrim(eventnum):end);
    AXID1.LAX=AXID1.LAX(ttrim(eventnum):end);
    AXID1.LAY=AXID1.LAY(ttrim(eventnum):end);
    AXID1.p_time=AXID1.p_time(ptrim(eventnum):end);
    AXID1.BDO=AXID1.BDO(ptrim(eventnum):end);
    
    AXID1.t_time(272983:273043)=[];AXID1.LAX(272983:273043)=[];AXID1.LAY(272983:273043)=[]; % jump here
    AXID1.LAX(272983:end)=AXID1.LAX(272983:end)-(AXID1.LAX(272983)-AXID1.LAX(272982));
    AXID1.LAY(272983:end)=AXID1.LAY(272983:end)-(AXID1.LAY(272983)-AXID1.LAY(272982));
    
    AXID1.LAX(259904:end)=AXID1.LAX(259904:end)-(AXID1.LAX(259904)-AXID1.LAX(259903));
    AXID1.LAY(259904:end)=AXID1.LAY(259904:end)-(AXID1.LAY(259904)-AXID1.LAY(259903));
    
    AXID1.t_time(172600:172700)=[];AXID1.LAX(172600:172700)=[];AXID1.LAY(172600:172700)=[]; % jump here
    AXID1.LAX(172600:end)=AXID1.LAX(172600:end)-(AXID1.LAX(172600)-AXID1.LAX(172599));
    AXID1.LAY(172600:end)=AXID1.LAY(172600:end)-(AXID1.LAY(172600)-AXID1.LAY(172599));
elseif eventnum==2 % no ASHES data
    AXCC1.t_time=AXCC1.t_time(ttrim(eventnum):end);
    AXCC1.LAX=AXCC1.LAX(ttrim(eventnum):end);
    AXCC1.LAY=AXCC1.LAY(ttrim(eventnum):end);
    AXCC1.p_time=AXCC1.p_time(ptrim(eventnum):end);
    AXCC1.BDO=AXCC1.BDO(ptrim(eventnum):end);
    
    AXEC2.t_time=AXEC2.t_time(ttrim(eventnum):end);
    AXEC2.LAX=AXEC2.LAX(ttrim(eventnum):end);
    AXEC2.LAY=AXEC2.LAY(ttrim(eventnum):end);
    AXEC2.p_time=AXEC2.p_time(ptrim(eventnum):end);
    AXEC2.BDO=AXEC2.BDO(ptrim(eventnum):end);
    
    AXID1.t_time=AXID1.time(ttrim(eventnum):end);
    AXID1.LAX=AXID1.LAX(ttrim(eventnum):end);
    AXID1.LAY=AXID1.LAY(ttrim(eventnum):end);
    AXID1.p_time=AXID1.p_time(ptrim(eventnum):end);
    AXID1.BDO=AXID1.BDO(ptrim(eventnum):end);
elseif eventnum==3 % full ASHES data
    AXCC1.t_time=AXCC1.t_time(ttrim(eventnum):end);
    AXCC1.LAX=AXCC1.LAX(ttrim(eventnum):end);
    AXCC1.LAY=AXCC1.LAY(ttrim(eventnum):end);
    AXCC1.p_time=AXCC1.p_time(ptrim(eventnum):end);
    AXCC1.BDO=AXCC1.BDO(ptrim(eventnum):end);
    
    AXCC1.t_time(401450:401510)=[];AXCC1.LAX(401450:401510)=[];AXCC1.LAY(401450:401510)=[]; % jump here
    AXCC1.LAX(401450:end)=AXCC1.LAX(401450:end)-(AXCC1.LAX(401450)-AXCC1.LAX(401449));
    AXCC1.LAY(401450:end)=AXCC1.LAY(401450:end)-(AXCC1.LAY(401450)-AXCC1.LAY(401449));
    
    AXCC1.t_time(369276:369932)=[];AXCC1.LAX(369276:369932)=[];AXCC1.LAY(369276:369932)=[]; % jump here
    AXCC1.LAX(369276:end)=AXCC1.LAX(369276:end)-(AXCC1.LAX(369276)-AXCC1.LAX(369275));
    AXCC1.LAY(369276:end)=AXCC1.LAY(369276:end)-(AXCC1.LAY(369276)-AXCC1.LAY(369275));
    
    AXCC1.t_time(352370:352430)=[];AXCC1.LAX(352370:352430)=[];AXCC1.LAY(352370:352430)=[]; % jump here
    AXCC1.LAX(352370:end)=AXCC1.LAX(352370:end)-(AXCC1.LAX(352370)-AXCC1.LAX(352369));
    AXCC1.LAY(352370:end)=AXCC1.LAY(352370:end)-(AXCC1.LAY(352370)-AXCC1.LAY(352369));
    
    AXCC1.t_time(328160:328610)=[];AXCC1.LAX(328160:328610)=[];AXCC1.LAY(328160:328610)=[]; % jump here
    AXCC1.LAX(328160:end)=AXCC1.LAX(328160:end)-(AXCC1.LAX(328160)-AXCC1.LAX(328159));
    AXCC1.LAY(328160:end)=AXCC1.LAY(328160:end)-(AXCC1.LAY(328160)-AXCC1.LAY(328159));
    
    AXCC1.t_time(170680:171525)=[];AXCC1.LAX(170680:171525)=[];AXCC1.LAY(170680:171525)=[]; % jump here
    AXCC1.LAX(170680:end)=AXCC1.LAX(170680:end)-(AXCC1.LAX(170680)-AXCC1.LAX(170679));
    AXCC1.LAY(170680:end)=AXCC1.LAY(170680:end)-(AXCC1.LAY(170680)-AXCC1.LAY(170679));
    
    AXEC2.LAY=AXEC2.LAY(159000:end);
    AXEC2.t_time=AXEC2.t_time(end-length(AXEC2.LAY)+1:end);
    AXEC2.LAX=AXEC2.LAX(end-length(AXEC2.LAY)+1:end);
    AXEC2.p_time=AXEC2.p_time(ptrim(eventnum):end);
    AXEC2.BDO=AXEC2.BDO(ptrim(eventnum):end);
    
    AXEC2.t_time(29150:29210)=[];AXEC2.LAX(29150:29210)=[];AXEC2.LAY(29150:29210)=[]; % jump here
    AXEC2.LAX(29150:end)=AXEC2.LAX(29150:end)-(AXEC2.LAX(29150)-AXEC2.LAX(29149));
    AXEC2.LAY(29150:end)=AXEC2.LAY(29150:end)-(AXEC2.LAY(29150)-AXEC2.LAY(29149));
    
    AXID1.t_time=AXID1.t_time(ttrim(eventnum):end);
    AXID1.LAX=AXID1.LAX(ttrim(eventnum):end);
    AXID1.LAY=AXID1.LAY(ttrim(eventnum):end);
    AXID1.p_time=AXID1.p_time(ptrim(eventnum):end);
    AXID1.BDO=AXID1.BDO(ptrim(eventnum):end);
    
    AXID1.t_time(400915:400975)=[];AXID1.LAX(400915:400975)=[];AXID1.LAY(400915:400975)=[]; % jump here
    AXID1.LAX(400915:end)=AXID1.LAX(400915:end)-(AXID1.LAX(400915)-AXID1.LAX(400914));
    AXID1.LAY(400915:end)=AXID1.LAY(400915:end)-(AXID1.LAY(400915)-AXID1.LAY(400914));
    
    AXID1.t_time(368470:368600)=[];AXID1.LAX(368470:368600)=[];AXID1.LAY(368470:368600)=[]; % jump here
    AXID1.LAX(368470:end)=AXID1.LAX(368470:end)-(AXID1.LAX(368470)-AXID1.LAX(368469));
    AXID1.LAY(368470:end)=AXID1.LAY(368470:end)-(AXID1.LAY(368470)-AXID1.LAY(368469));
    
    AXID1.t_time(162150:162210)=[];AXID1.LAX(162150:162210)=[];AXID1.LAY(162150:162210)=[]; % jump here
    AXID1.LAX(162150:end)=AXID1.LAX(162150:end)-(AXID1.LAX(162150)-AXID1.LAX(162149));
    AXID1.LAY(162150:end)=AXID1.LAY(162150:end)-(AXID1.LAY(162150)-AXID1.LAY(162149));
    
    AXID1.t_time(82830:82890)=[];AXID1.LAX(82830:82890)=[];AXID1.LAY(82830:82890)=[]; % jump here
    AXID1.LAX(82830:end)=AXID1.LAX(82830:end)-(AXID1.LAX(82830)-AXID1.LAX(82829));
    AXID1.LAY(82830:end)=AXID1.LAY(82830:end)-(AXID1.LAY(82830)-AXID1.LAY(82829));
    
    % ASHES
    ASHES.t_time=ASHES.t_time(ttrim(eventnum):end);
    ASHES.LAX=ASHES.LAX(ttrim(eventnum):end);
    ASHES.LAY=ASHES.LAY(ttrim(eventnum):end);
    ASHES.p_time=ASHES.p_time(ptrim(eventnum):end);
    ASHES.BDO=ASHES.BDO(ptrim(eventnum):end);
    
    ASHES.t_time(92315:92375)=[];ASHES.LAX(92315:92375)=[];ASHES.LAY(92315:92375)=[]; % jump here
    ASHES.LAX(92315:end)=ASHES.LAX(92315:end)-(ASHES.LAX(92315)-ASHES.LAX(92314));
    ASHES.LAY(92315:end)=ASHES.LAY(92315:end)-(ASHES.LAY(92315)-ASHES.LAY(92314));
    
    ASHES.t_time(25730:25930)=[];ASHES.LAX(25730:25930)=[];ASHES.LAY(25730:25930)=[]; % jump here
    ASHES.LAX(25730:end)=ASHES.LAX(25730:end)-(ASHES.LAX(25730)-ASHES.LAX(25729));
    ASHES.LAY(25730:end)=ASHES.LAY(25730:end)-(ASHES.LAY(25730)-ASHES.LAY(25729));
elseif eventnum==4 % full ASHES data
    % AXCC1
    AXCC1.t_time=AXCC1.t_time(ttrim(eventnum):end);
    AXCC1.LAX=AXCC1.LAX(ttrim(eventnum):end);
    AXCC1.LAY=AXCC1.LAY(ttrim(eventnum):end);
    AXCC1.p_time=AXCC1.p_time(ptrim(eventnum):end);
    AXCC1.BDO=AXCC1.BDO(ptrim(eventnum):end);
    
    % AXEC2
    AXEC2.t_time=AXEC2.t_time(ttrim(eventnum):end);
    AXEC2.LAX=AXEC2.LAX(ttrim(eventnum):end);
    AXEC2.LAY=AXEC2.LAY(ttrim(eventnum):end);
    AXEC2.p_time=AXEC2.p_time(ptrim(eventnum):end);
    AXEC2.BDO=AXEC2.BDO(ptrim(eventnum):end);
    
    AXEC2.t_time(170120:170180)=[];AXEC2.LAX(170120:170180)=[];AXEC2.LAY(170120:170180)=[]; % jump here
    AXEC2.LAX(170120:end)=AXEC2.LAX(170120:end)-(AXEC2.LAX(170120)-AXEC2.LAX(170119));
    AXEC2.LAY(170120:end)=AXEC2.LAY(170120:end)-(AXEC2.LAY(170120)-AXEC2.LAY(170119));
    
    AXEC2.t_time(98930:98990)=[];AXEC2.LAX(98930:98990)=[];AXEC2.LAY(98930:98990)=[]; % jump here
    AXEC2.LAX(98930:end)=AXEC2.LAX(98930:end)-(AXEC2.LAX(98930)-AXEC2.LAX(98929));
    AXEC2.LAY(98930:end)=AXEC2.LAY(98930:end)-(AXEC2.LAY(98930)-AXEC2.LAY(98929));
    
    % AXID1
    AXID1.t_time=AXID1.t_time(ttrim(eventnum):end);
    AXID1.LAX=AXID1.LAX(ttrim(eventnum):end);
    AXID1.LAY=AXID1.LAY(ttrim(eventnum):end);
    AXID1.p_time=AXID1.p_time(ptrim(eventnum):end);
    AXID1.BDO=AXID1.BDO(ptrim(eventnum):end);
    
    AXID1.t_time(35670:35730)=[];AXID1.LAX(35670:35730)=[];AXID1.LAY(35670:35730)=[]; % jump here
    AXID1.LAX(35670:end)=AXID1.LAX(35670:end)-(AXID1.LAX(35670)-AXID1.LAX(35669));
    AXID1.LAY(35670:end)=AXID1.LAY(35670:end)-(AXID1.LAY(35670)-AXID1.LAY(35669));
    
    % ASHES
    ASHES.t_time=ASHES.t_time(ttrim(eventnum):end);
    ASHES.LAX=ASHES.LAX(ttrim(eventnum):end);
    ASHES.LAY=ASHES.LAY(ttrim(eventnum):end);
    ASHES.p_time=ASHES.p_time(ptrim(eventnum):end);
    ASHES.BDO=ASHES.BDO(ptrim(eventnum):end);
    
    ASHES.t_time(39280:39340)=[];ASHES.LAX(39280:39340)=[];ASHES.LAY(39280:39340)=[]; % jump here
    ASHES.LAX(39280:end)=ASHES.LAX(39280:end)-(ASHES.LAX(39280)-ASHES.LAX(39279));
    ASHES.LAY(39280:end)=ASHES.LAY(39280:end)-(ASHES.LAY(39280)-ASHES.LAY(39279));
elseif eventnum==5 % full ASHES data
    % AXCC1
    AXCC1.t_time=AXCC1.t_time(ttrim(eventnum):end);
    AXCC1.LAX=AXCC1.LAX(ttrim(eventnum):end);
    AXCC1.LAY=AXCC1.LAY(ttrim(eventnum):end);
    AXCC1.p_time=AXCC1.p_time(ptrim(eventnum):end);
    AXCC1.BDO=AXCC1.BDO(ptrim(eventnum):end);
    
    AXCC1.t_time(74600:74660)=[];AXCC1.LAX(74600:74660)=[];AXCC1.LAY(74600:74660)=[]; % jump here
    AXCC1.LAX(74600:end)=AXCC1.LAX(74600:end)-(AXCC1.LAX(74600)-AXCC1.LAX(74599));
    AXCC1.LAY(74600:end)=AXCC1.LAY(74600:end)-(AXCC1.LAY(74600)-AXCC1.LAY(74599));
    
    % AXEC2
    AXEC2.t_time=AXEC2.t_time(ttrim(eventnum):end);
    AXEC2.LAX=AXEC2.LAX(ttrim(eventnum):end);
    AXEC2.LAY=AXEC2.LAY(ttrim(eventnum):end);
    AXEC2.p_time=AXEC2.p_time(ptrim(eventnum):end);
    AXEC2.BDO=AXEC2.BDO(ptrim(eventnum):end);
    
    AXEC2.t_time(129530:129590)=[];AXEC2.LAX(129530:129590)=[];AXEC2.LAY(129530:129590)=[]; % jump here
    AXEC2.LAX(129530:end)=AXEC2.LAX(129530:end)-(AXEC2.LAX(129530)-AXEC2.LAX(129529));
    AXEC2.LAY(129530:end)=AXEC2.LAY(129530:end)-(AXEC2.LAY(129530)-AXEC2.LAY(129529));
    
    AXEC2.t_time(104250:104310)=[];AXEC2.LAX(104250:104310)=[];AXEC2.LAY(104250:104310)=[]; % jump here
    AXEC2.LAX(104250:end)=AXEC2.LAX(104250:end)-(AXEC2.LAX(104250)-AXEC2.LAX(104249));
    AXEC2.LAY(104250:end)=AXEC2.LAY(104250:end)-(AXEC2.LAY(104250)-AXEC2.LAY(104249));
    
    % AXID1
    AXID1.t_time=AXID1.t_time(ttrim(eventnum):end);
    AXID1.LAX=AXID1.LAX(ttrim(eventnum):end);
    AXID1.LAY=AXID1.LAY(ttrim(eventnum):end);
    AXID1.p_time=AXID1.p_time(ptrim(eventnum):end);
    AXID1.BDO=AXID1.BDO(ptrim(eventnum):end);
    
    AXID1.t_time(61190:61250)=[];AXID1.LAX(61190:61250)=[];AXID1.LAY(61190:61250)=[]; % jump here
    AXID1.LAX(61190:end)=AXID1.LAX(61190:end)-(AXID1.LAX(61190)-AXID1.LAX(61189));
    AXID1.LAY(61190:end)=AXID1.LAY(61190:end)-(AXID1.LAY(61190)-AXID1.LAY(61189));
    
    AXID1.t_time(28680:28740)=[];AXID1.LAX(28680:28740)=[];AXID1.LAY(28680:28740)=[]; % jump here
    AXID1.LAX(28680:end)=AXID1.LAX(28680:end)-(AXID1.LAX(28680)-AXID1.LAX(28679));
    AXID1.LAY(28680:end)=AXID1.LAY(28680:end)-(AXID1.LAY(28680)-AXID1.LAY(28679));
    
    % ASHES
    ASHES.t_time=ASHES.t_time(ttrim(eventnum):end);
    ASHES.LAX=ASHES.LAX(ttrim(eventnum):end);
    ASHES.LAY=ASHES.LAY(ttrim(eventnum):end);
    ASHES.p_time=ASHES.p_time(ptrim(eventnum):end);
    ASHES.BDO=ASHES.BDO(ptrim(eventnum):end);
    
    ASHES.t_time(197920:197980)=[];ASHES.LAX(197920:197980)=[];ASHES.LAY(197920:197980)=[]; % jump here
    ASHES.LAX(197920:end)=ASHES.LAX(197920:end)-(ASHES.LAX(197920)-ASHES.LAX(197919));
    ASHES.LAY(197920:end)=ASHES.LAY(197920:end)-(ASHES.LAY(197920)-ASHES.LAY(197919));
end

save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXEC2clean'],'AXEC2')
save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXCC1clean'],'AXCC1')
save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXID1clean'],'AXID1')
if ~exclude_ashes
    save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/ASHESclean'],'ASHES')
end