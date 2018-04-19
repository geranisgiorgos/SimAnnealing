function [minimum,fval] = sm_anneal(p, parent, A, C, eqin, foit, kath, options)
%  Είσοδοι:
%  P: η απολαβή των φοιτητών αν πάρουν πτυχιακή σε ένα συγκεκριμένο
%  καθηγητή
%
%  PARENT: η αρχική εκτίμηση
% 
%  A: ο πίνακας των περιορισμών
% 
%  C: ο πίνακας με τις τιμές των περιορισμών
% 
%  EQIN: ο πίνακας με το είδος των περιορισμών
% 
%  FOIT: ο αριθμός των φοιτητών
% 
%  KATH: ο αριθμός των καθηγητών
%
%  OPTIONS: η δομή με τις επιλογές για τον αλγόριθμο. Οι επιλογές της είναι
%  οι εξής:
%       Verbosity: Controls output to the screen.
%                  0 suppresses all output
%                  1 gives final report only [default]
%                  2 gives temperature changes and final report 
%       Generator: Generates a new solution from an old one.
%                  Any function handle that takes a solution as input and
%                  gives a valid solution (i.e. some point in the solution
%                  space) as output.
%                  The default function generates a row vector which
%                  slightly differs from the input vector in one element:
%                  @(x) (x+(randperm(length(x))==length(x))*randn/100)
%                  Other examples of possible solution generators:
%                  @(x) (rand(3,1)): Picks a random point in the unit cube
%                  @(x) (ceil([9 5].*rand(2,1))): Picks a point in a 9-by-5
%                                                 discrete grid
%                  Note that if you use the default generator, ANNEAL only
%                  works on row vectors. For loss functions that operate on
%                  column vectors, use this generator instead of the
%                  default:
%                  @(x) (x(:)'+(randperm(length(x))==length(x))*randn/100)'
%        InitTemp: The initial temperature, can be any positive number.
%                  Default is 1.
%        StopTemp: Temperature at which to stop, can be any positive number
%                  smaller than InitTemp. 
%                  Default is 1e-8.
%         StopVal: Value at which to stop immediately, can be any output of
%                  LOSS that is sufficiently low for you.
%                  Default is -Inf.
%       CoolSched: Generates a new temperature from the previous one.
%                  Any function handle that takes a scalar as input and
%                  returns a smaller but positive scalar as output. 
%                  Default is @(T) (.8*T)
%      MaxConsRej: Maximum number of consecutive rejections, can be any
%                  positive number.
%                  Default is 1000.
%        MaxTries: Maximum number of tries within one temperature, can be
%                  any positive number.
%                  Default is 300.
%      MaxSuccess: Maximum number of successes within one temperature, can
%                  be any positive number.
%                  Default is 20.

def = struct(...
        'CoolSched',@(T) (.8*T),...
        'Generator',@(x) (x+(randperm(length(x))==length(x))*randn/100),...
        'InitTemp',1,...
        'MaxConsRej',1000,...
        'MaxSuccess',20,...
        'MaxTries',300,...
        'StopTemp',1e-8,...
        'StopVal',+Inf,...
        'Verbosity',1);

% Check input
if nargin<7, %user input error
    error('MATLAB:anneal:Syntax error','You need to input at least the five parameters.');
elseif nargin<8, %user gave no options structure, use default
    options=def;
else %user gave all input, check if options structure is complete
    if ~isstruct(options)
        error('MATLAB:anneal:badOptions',...
            'Input argument ''options'' is not a structure')
    end
    fs = {'CoolSched','Generator','InitTemp','MaxConsRej',...
        'MaxSuccess','MaxTries','StopTemp','StopVal','Verbosity'};
    for nm=1:length(fs)
        if ~isfield(options,fs{nm}), options.(fs{nm}) = def.(fs{nm}); end
    end
end

% main settings
newsol = options.Generator;      % neighborhood space function
Tinit = options.InitTemp;        % initial temp
minT = options.StopTemp;         % stopping temp
cool = options.CoolSched;        % annealing schedule
maxF = options.StopVal;
max_consec_rejections = options.MaxConsRej;
max_try = options.MaxTries;
max_success = options.MaxSuccess;
report = options.Verbosity;
k = 1;                           % boltzmann constant

% counters etc
itry = 0;
success = 0;
finished = 0;
consec = 0;
T = Tinit;
initenergy = sum(sum(parent.*p));
oldenergy = initenergy;
total = 0;
if report==2, fprintf(1,'\n  T = %7.5f, loss = %10.5f\n',T,oldenergy); end
while ~finished
    itry = itry+1; % just an iteration counter
    current = parent; 
    
    % % Stop / decrement T criteria
    if itry >= max_try || success >= max_success;
        if T > minT || consec >= max_consec_rejections;
            finished = 1;
            total = total + itry;
            break;
        else
            T = cool(T);  % decrease T according to cooling schedule
            if report==2, % output
                fprintf(1,'  T = %7.5f, loss = %10.5f\n',T,oldenergy);
            end
            total = total + itry;
            itry = 1;
            success = 1;
        end
    end
    
    mesos=sum(C(foit+1:length(C)))/foit;
    anatfoitites = zeros(1,foit);
    newparam = zeros(1, foit*kath);
    foitites = [1:1:foit];
    synoloptyx = foit;
    newparam = zeros(1, foit*kath);
    A_new = zeros(foit, kath);
    
    % παραγωγή νέας τυχαίας θερμοκρασίας
    kathigitis = 1;
    while(synoloptyx > 0)
        temp = ceil(rand(1)*length(foitites));
        newparam((foitites(temp)-1)*kath + kathigitis) = 1;
        A_new(foitites(temp), kathigitis) = 1;
        foitites = [foitites(1:temp-1) foitites(temp+1:length(foitites))];
        kathigitis = kathigitis + 1;
        if(kathigitis > kath)
            kathigitis = 1;
        end
        synoloptyx = synoloptyx - 1;
    end
    newenergy = sum(sum(A_new.*p));
    
    %έλεγχος αν η νέα λύση ικανοποιεί τους περιορισμούς
    flag = 0;
    temp = A * newparam';  %πολλαπλασιασμός του A*x
    for ii=1:length(temp)
        if(eqin(ii) == 0)  %έλεγχος για τους ισοτικούς περιορισμούς
            if(abs(temp(ii) - C(ii)) > 1e-6)
                flag = 1;
                break;
            end
        else if(eqin(ii) == -1) %έλεγχους για τους ανισοτικούς περιορισμούς τύπου <
                if(temp(ii) - C(ii) > 1e-6)
                    flag = 1;
                    break;
                end
            else if(eqin(ii) == 1) %έλεγχους για τους ανισοτικούς περιορισμούς τύπου >
                    if(temp(ii) - C(ii) < 1e-6)
                        flag = 1;
                        break;
                    end   
                end
            end
        end
    end
    %αν δεν ικανοποιείται κάποιος περιορισμός τότε δεν κρατάμε τη λύση και
    %προχωράμε στην επόμενη επανάληψη
    if(flag == 1)
        continue;
    end
    
    if (newenergy > maxF),
        parent = newparam; 
        oldenergy = newenergy;
        break
    end
    
    if (oldenergy-newenergy < 1e-6)
        parent = newparam;
        oldenergy = newenergy;
        success = success+1;
        consec = 0;
    else
        if (rand > exp( (oldenergy-newenergy)/(k*T) ));
            parent = newparam;
            oldenergy = newenergy;
            success = success+1;
        else
            consec = consec+1;
        end
    end
end

minimum = parent;
fval = oldenergy;

%if report;
%    fprintf(1, '\n  Initial temperature:     \t%g\n', Tinit);
%    fprintf(1, '  Final temperature:       \t%g\n', T);
%    fprintf(1, '  Consecutive rejections:  \t%i\n', consec);
%    fprintf(1, '  Number of function calls:\t%i\n', total);
%    fprintf(1, '  Total final loss:        \t%g\n', fval);
%end

