function run_chain(i,x0,sv_file,get_err,proprnd,thin,nsamples,nsegments)
    % Function that runs the MH chains and saves the files

    % i - chain number
    % x0 - intial parameter guess
    % sv_file - the file name to save this chain too
    % get_err - get LL function
    % proprnd - proposal distribution function for new parameter guesses
    % thin - thinning rate, only one in this value will be saved
    % nsamples - number of samples per each chain
    % nsegments - number of segments to break each chain into

    pause(i)
    rng('shuffle')
    for iq = 1:nsegments
        
        [mh_smpli,accepti,mh_valuei] = MetHast(x0,nsamples,'logpdf',get_err,'proprnd', ...
            proprnd,'symmetric',1,'thin',thin);
        
        eval(['mh_smpl_',num2str(iq),'= mh_smpli(1:end-1,:);']);
        eval(['mh_value_',num2str(iq),'= mh_valuei(1:end-1,:);']);
        
        x0 = mh_smpli(end,:);
        
        accepti
        if iq==1
            save(sv_file,['mh_value_',num2str(iq)],['mh_smpl_',num2str(iq)]);
        else
            save(sv_file,['mh_value_',num2str(iq)],['mh_smpl_',num2str(iq)],'-append');
        end
    end
end
