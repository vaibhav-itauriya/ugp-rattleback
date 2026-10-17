% =========================================================================
% RATTLEBACK ON A ROUGH PLANE - SLIPPING CONTACT WITH FRICTION
% =========================================================================
% Companion code to: V. Itauriya and D. Mukherjee, "Slipping-contact friction
% in rattleback dynamics: extension of a rolling formulation using the
% Moffatt-Shimomura contact laws".
% Rolling formulation: D. Mukherjee, "On preferential rotation of the
%   rattleback" (1-2-3 Euler angles, rolling without slipping, 6 ODEs).
% Tangential contact laws: H. K. Moffatt and Y. Shimomura, Nature 416, 385
%   (2002).
%
% Requires MATLAB R2021a or later (name=value args, exportgraphics, tiles).
%
% -------------------------------------------------------------------------
% WHAT THIS FILE RUNS   (switch studies on/off in config() below)
% -------------------------------------------------------------------------
%   MAIN       four contact models x two senses of initial spin
%              -> per-model omega / theta / contact figures
%              -> two cross-model comparison figures
%   EPS SWEEP  regularisation convergence: does the reversal time settle as
%              eps -> 0?  This is the largest stated uncertainty.
%   RHO SWEEP  contact patch radius: locates the physically realistic range
%              and shows where drilling friction kills the effect entirely.
%   MU SWEEP   rolling-limit convergence: Coulomb must reproduce the rolling
%              benchmark as mu grows. Strongest single validation.
%   W0 COMPARE two initial spins, to check the dependence of the reversal
%              time on omega_0.
%   FIXED N    repeats the Coulomb runs with N held at m*g, to measure the
%              error of that approximation against the solved N.
%   ANIMATE    writes one GIF per contact model.
%
% Figures are exported to ./figures at 300 dpi when CFG.export is true.
%
% -------------------------------------------------------------------------
% FORMULATION
% -------------------------------------------------------------------------
% Friction is a function of the SLIP velocity v_P. Under rolling v_P = 0, so
% friction cannot be added to the rolling model - the constraint has to be
% released first. That makes v_G independent (3 extra states, 9 in total)
% and means the contact force can no longer be eliminated. Its normal part
% N is NOT m*g: in the rolling benchmark N runs from about 0.20*m*g to
% 5.9*m*g. N is solved from the contact-maintenance condition:
%
%     zG      = h(n) = sqrt(sum a_i^2 n_i^2) - b*n3 = -r.n
%     zGddot  = [-rdot.ndot + r.(w x ndot)] + wdot.(n x r) = c + d.wdot
%     N - m*g = m*zGddot                     (friction is tangential)
%     I*wdot  = r x (N*n + F) - w x I*w
%  => wdot = A0 + N*(B + A1),  N = m*(g + c + d.A0)/(1 - m*d.B - m*d.A1)
%
% one scalar division per RHS evaluation, no iteration, no assumption.
% =========================================================================

function z = rattleback_simulator_123H()
clear all; close all; clc; z = 0;

CFG = config();
P0  = parameters();
MDL = {};

if CFG.export, CFG.figdir = prepare_figdir(CFG.figdir); end

if CFG.main,      MDL = study_main(P0,CFG); end
if CFG.epsSweep,  study_eps(P0,CFG);        end
if CFG.rhoSweep,  study_rho(P0,CFG);        end
if CFG.muSweep,   study_mu(P0,CFG);         end
if CFG.w0Compare, study_w0(P0,CFG);         end
if CFG.fixNCompare, study_fixN(P0,CFG);     end

if CFG.animate && ~isempty(MDL)
    for k = 1:numel(MDL)
        animate(MDL{k}.T, MDL{k}.Y2, MDL{k}.P, MDL{k}.tag, CFG);
    end
end
fprintf('\nDone.\n');
end

% =========================================================================
% CONFIGURATION
% =========================================================================
function CFG = config()
CFG.main      = true;
CFG.epsSweep  = true;
CFG.rhoSweep  = true;
CFG.muSweep   = true;
CFG.w0Compare = true;
CFG.fixNCompare = true;     % solved N against N = m*g
CFG.animate   = false;      % slow: roughly a minute per model

CFG.w0     = 5.0;           % initial spin magnitude              [rad/s]
CFG.tEnd   = 8.0;           % simulated time                          [s]
CFG.nOut   = 20001;         % output points (does NOT set accuracy)
CFG.nSweep = 4001;          % lighter grid for the parameter sweeps
CFG.relTol = 1e-8;
CFG.absTol = 1e-10;

CFG.export     = true;
CFG.figdir     = 'figures';   % created next to this .m file
CFG.dpi        = 300;
CFG.animFrames = 200;

% ---- house style --------------------------------------------------------
CFG.fs   = 10;
CFG.font = 'Helvetica';

% Okabe-Ito derived palette: separable in greyscale and for colour-blind
% readers.
CFG.col.ink     = [0.12 0.15 0.19];   % rolling / neutral
CFG.col.blue    = [0.00 0.45 0.70];   % Coulomb        / omega_1
CFG.col.verm    = [0.84 0.37 0.00];   % Coulomb + spin / omega_2
CFG.col.green   = [0.00 0.62 0.45];   % omega_3
CFG.col.purple  = [0.49 0.18 0.56];   % viscous
CFG.col.grey    = [0.55 0.60 0.65];
CFG.col.pref    = [0.00 0.45 0.70];
CFG.col.nonpref = [0.84 0.37 0.00];
end

% =========================================================================
% PHYSICAL PARAMETERS (Kane and Levinson, Int. J. Non-Linear Mech. 17, 175, 1982)
% =========================================================================
function P = parameters()
P.m  = 1;                                   % Mass                    [kg]
P.g  = 9.81;                                % Gravity              [m/s^2]
P.a1 = 0.2; P.a2 = 0.03; P.a3 = 0.02;       % Semi-axes                [m]
P.b  = 0.005;                               % CM depth below centre    [m]
Ixx  = 0.0002; Iyy = 0.0016; Izz = 0.00175;
Ixy  = -0.00005;    % DO NOT SET TO ZERO: the chirality, and with it the
                    % reversal, dies along with this one term
P.I  = [Ixx Ixy 0; Ixy Iyy 0; 0 0 Izz];
P.a  = [P.a1; P.a2; P.a3];

P.law = 'coulomb';
P.mu  = 0.15;       % Coulomb coefficient                              [-]
P.eps = 3e-5;       % Slip regularisation                            [m/s]
                    % Results are converged in eps at this value: the
                    % reversal time changes by under 1% between 3e-5 and 1e-5
                    % (see study_eps). Values near 1e-3 are NOT converged.
P.rho = 0.0;        % Contact patch radius (drilling friction)         [m]
P.cv  = 500;        % Viscous coefficient                           [kg/s]
P.kb  = 0;          % Baumgarte gain, 0 = off                        [1/s]
P.fixN = false;     % true holds N at m*g (comparison only - this violates
                    % the contact condition; see study_fixN)
end

% =========================================================================
% STUDY 1 - MAIN
% =========================================================================
function MDL = study_main(P0,CFG)
banner('MAIN STUDY - four contact models, both senses of initial spin');
t = linspace(0,CFG.tEnd,CFG.nOut);

MDL = {};
Pr = P0; Pr.law = 'rolling';
MDL{end+1} = struct('P',Pr,'name','Rolling', ...
                    'tag','rolling','col',CFG.col.ink);

Pc = P0; Pc.law = 'coulomb'; Pc.rho = 0;
MDL{end+1} = struct('P',Pc,'name',sprintf('Coulomb, \\mu = %.2f',Pc.mu), ...
                    'tag','coulomb','col',CFG.col.blue);

% Effective drilling lever arm rho = (3*pi/16)*c for a fully sliding Hertzian
% patch of radius c. With radii of curvature a1^2/a3 = 2.0 m and
% a2^2/a3 = 0.045 m at the contact, c is ~0.27 mm for steel on steel and
% ~0.88 mm for a polymer body on glass, giving rho ~0.16-0.52 mm. The upper
% value is used here. At rho = 2 mm the torque stops the body before the
% instability can develop (see study_rho).
Ps = P0; Ps.law = 'coulomb'; Ps.rho = 0.5e-3;
MDL{end+1} = struct('P',Ps,'name',sprintf('Coulomb + spin, \\rho = %.1f mm',Ps.rho*1e3), ...
                    'tag','coulomb_spin','col',CFG.col.verm);

Pv = P0; Pv.law = 'viscous';
MDL{end+1} = struct('P',Pv,'name',sprintf('Viscous, c = %g',Pv.cv), ...
                    'tag','viscous','col',CFG.col.purple);

hdr_table();
for k = 1:numel(MDL)
    Pk = MDL{k}.P;
    [T,Y1,D1] = run_case(t,Pk,+CFG.w0,CFG);
    [~,Y2,D2] = run_case(t,Pk,-CFG.w0,CFG);
    MDL{k}.T = T; MDL{k}.Y1 = Y1; MDL{k}.Y2 = Y2;
    MDL{k}.D1 = D1; MDL{k}.D2 = D2;

    row_table(MDL{k}.name,'preferred',    T,Y1,D1,+CFG.w0);
    row_table(MDL{k}.name,'non-preferred',T,Y2,D2,-CFG.w0);

    fig_omega(T,Y1,Y2,MDL{k},CFG);
    fig_theta(T,Y1,Y2,MDL{k},CFG);
    fig_contact(D1,D2,MDL{k},CFG);
end

fig_compare_contact(MDL,CFG);
fig_compare_omega(MDL,CFG);
end

% =========================================================================
% STUDY 2 - REGULARISATION CONVERGENCE
% The reversal time depends on eps. A converged value has to be quoted, not
% a single run. This is the largest single uncertainty in the results.
% =========================================================================
function study_eps(P0,CFG)
banner('EPSILON SWEEP - regularisation convergence');
epsList = [1e-2 3e-3 1e-3 3e-4 1e-4 3e-5 1e-5];
t  = linspace(0,CFG.tEnd,CFG.nSweep);
tr = nan(size(epsList)); we = nan(size(epsList));

fprintf('%12s %12s %12s %12s\n','eps (m/s)','t_rev (s)','w3(end)','max|gv|');
fprintf('%s\n',repmat('-',1,52));
for k = 1:numel(epsList)
    P = P0; P.law = 'coulomb'; P.rho = 0; P.eps = epsList(k);
    [T,Y,D] = run_case(t,P,-CFG.w0,CFG);
    R = reversal(T,Y(:,3),-CFG.w0);
    tr(k) = R.t_cross; we(k) = Y(end,3);
    fprintf('%12.1e %12s %12.4f %12.2e\n', ...
            epsList(k), fmt_t(R.t_cross), we(k), max(abs(D.gv)));
end

f  = newfig('eps convergence',[620 420],CFG);
tl = tiledlayout(f,2,1,'Padding','compact','TileSpacing','compact');

ax = nexttile(tl);
semilogx(ax,epsList,tr,'o-','Color',CFG.col.blue,'LineWidth',1.6, ...
         'MarkerFaceColor','w','MarkerSize',6);
style(ax,CFG);
xlabel(ax,'$\varepsilon$ (m/s)','Interpreter','latex');
ylabel(ax,'$t_{\mathrm{rev}}$ (s)','Interpreter','latex');
title(ax,'Reversal time vs slip regularisation');

ax = nexttile(tl);
semilogx(ax,epsList,we,'s-','Color',CFG.col.verm,'LineWidth',1.6, ...
         'MarkerFaceColor','w','MarkerSize',6);
style(ax,CFG);
xlabel(ax,'$\varepsilon$ (m/s)','Interpreter','latex');
ylabel(ax,'$\omega_3(t_{\mathrm{end}})$ (rad/s)','Interpreter','latex');
title(ax,'Recovered spin vs slip regularisation');

savefig_(f,'sweep_epsilon',CFG);
fprintf(['\nThe reversal time is converged where the upper curve is flat\n' ...
         '(eps <= 3e-5 m/s for these parameters).\n']);
end

% =========================================================================
% STUDY 3 - CONTACT PATCH RADIUS
% =========================================================================
function study_rho(P0,CFG)
banner('RHO SWEEP - drilling friction vs contact patch radius');
rhoList = [0 0.1 0.2 0.5 1.0 2.0]*1e-3;
t = linspace(0,CFG.tEnd,CFG.nSweep);
w3p = nan(size(rhoList)); w3n = nan(size(rhoList)); thn = nan(size(rhoList));

fprintf('%10s %12s %16s %18s %16s\n','rho (mm)','t_cross (s)', ...
        'w3 end (pref)','w3 end (non-pref)','theta max (deg)');
fprintf('%s\n',repmat('-',1,76));
for k = 1:numel(rhoList)
    P = P0; P.law = 'coulomb'; P.rho = rhoList(k);
    [~,Yp,~] = run_case(t,P,+CFG.w0,CFG);
    [Tn,Yn,~] = run_case(t,P,-CFG.w0,CFG);
    w3p(k) = Yp(end,3); w3n(k) = Yn(end,3);
    thn(k) = max(abs(Yn(:,5)))*180/pi;
    R = reversal(Tn,Yn(:,3),-CFG.w0);
    ts = fmt_t(R.t_cross); if ~isnan(R.t_cross) && ~R.true, ts = [ts '*']; end
    fprintf('%10.2f %12s %16.3f %18.3f %16.3f\n', ...
            rhoList(k)*1e3, ts, w3p(k), w3n(k), thn(k));
end

f  = newfig('rho sweep',[620 420],CFG);
tl = tiledlayout(f,2,1,'Padding','compact','TileSpacing','compact');

ax = nexttile(tl); hold(ax,'on');
hA = plot(ax,rhoList*1e3,w3p,'o-','Color',CFG.col.blue, ...
          'LineWidth',1.6,'MarkerFaceColor','w','MarkerSize',6);
hB = plot(ax,rhoList*1e3,w3n,'s-','Color',CFG.col.verm, ...
          'LineWidth',1.6,'MarkerFaceColor','w','MarkerSize',6);
hold(ax,'off'); style(ax,CFG);
legend(ax,[hA hB],{'preferred','non-preferred'},'Location','northeast');
xlabel(ax,'$\rho$ (mm)','Interpreter','latex');
ylabel(ax,'$\omega_3(t_{\mathrm{end}})$ (rad/s)','Interpreter','latex');
title(ax,'Spin-down vs contact patch radius');

ax = nexttile(tl);
plot(ax,rhoList*1e3,thn,'^-','Color',CFG.col.ink, ...
     'LineWidth',1.6,'MarkerFaceColor','w','MarkerSize',6);
style(ax,CFG);
xlabel(ax,'$\rho$ (mm)','Interpreter','latex');
ylabel(ax,'$\theta_{\max}$ (deg)','Interpreter','latex');
title(ax,'Wobble amplitude: the effect is extinguished near \rho = 2 mm');

savefig_(f,'sweep_rho',CFG);
fprintf(['\nHertzian estimate of the lever arm: rho = 0.16-0.52 mm (steel on\n' ...
         'steel to polymer on glass). At 2 mm the drilling torque stops the\n' ...
         'body before the instability develops. * = spin-down, not reversal.\n']);
end

% =========================================================================
% STUDY 4 - ROLLING-LIMIT CONVERGENCE
% As mu grows the Coulomb solution must reproduce the rolling benchmark.
% This exercises the geometry, the normal-force solve and the friction
% implementation at once, so it is the strongest single validation.
% =========================================================================
function study_mu(P0,CFG)
banner('MU SWEEP - convergence to the rolling benchmark');
muList = [0.05 0.15 0.5 1.0 2.0];
t = linspace(0,CFG.tEnd,CFG.nSweep);

Pr = P0; Pr.law = 'rolling';
[Tr,Yr,Dr] = run_case(t,Pr,-CFG.w0,CFG);
Rr = reversal(Tr,Yr(:,3),-CFG.w0);
fprintf('rolling benchmark: t_rev = %s s, N in [%.2f, %.2f] N\n\n', ...
        fmt_t(Rr.t_cross), min(Dr.N), max(Dr.N));

fprintf('%8s %12s %12s %18s\n','mu','t_rev (s)','w3(end)','N range (N)');
fprintf('%s\n',repmat('-',1,54));

f  = newfig('mu sweep',[680 400],CFG);
ax = axes('Parent',f); hold(ax,'on');
hh = gobjects(numel(muList)+1,1);
nm = cell(numel(muList)+1,1);
hh(1) = plot(ax,Tr,Yr(:,3),'-','Color',CFG.col.ink,'LineWidth',2.2);
nm{1} = 'Rolling benchmark';
sh = linspace(0.25,1.0,numel(muList));
for k = 1:numel(muList)
    P = P0; P.law = 'coulomb'; P.rho = 0; P.mu = muList(k);
    [T,Y,D] = run_case(t,P,-CFG.w0,CFG);
    R = reversal(T,Y(:,3),-CFG.w0);
    fprintf('%8.2f %12s %12.3f   [%6.2f, %7.2f]\n', ...
            muList(k), fmt_t(R.t_cross), Y(end,3), min(D.N), max(D.N));
    cc = CFG.col.blue*sh(k) + (1-sh(k))*[0.80 0.86 0.92];
    hh(k+1) = plot(ax,T,Y(:,3),'-','Color',cc,'LineWidth',1.3);
    nm{k+1} = sprintf('\\mu = %.2f',muList(k));
end
plot(ax,[t(1) t(end)],[0 0],':','Color',CFG.col.grey, ...
     'LineWidth',0.8,'HandleVisibility','off');
hold(ax,'off'); style(ax,CFG);
legend(ax,hh,nm,'Location','southeast');           % explicit handles
xlabel(ax,'$t$ (s)','Interpreter','latex');
ylabel(ax,'$\omega_3$ (rad/s)','Interpreter','latex');
title(ax,'Coulomb friction converges to the rolling solution as \mu grows');
savefig_(f,'sweep_mu',CFG);
end

% =========================================================================
% STUDY 5 - DEPENDENCE ON INITIAL SPIN
% In the rolling formulation a smaller initial spin in the unstable sense
% reverses later, because the rocking instability grows more slowly. Check
% whether Coulomb friction preserves that ordering.
% =========================================================================
function study_w0(P0,CFG)
banner('INITIAL SPIN COMPARISON');
w0List = [2.5 5.0];
t = linspace(0,CFG.tEnd,CFG.nSweep);

f  = newfig('initial spin',[680 440],CFG);
tl = tiledlayout(f,2,1,'Padding','compact','TileSpacing','compact');
names = {'Rolling','Coulomb, \mu = 0.15'};
laws  = {'rolling','coulomb'};
cols  = {CFG.col.ink, CFG.col.blue};

fprintf('%24s %8s %12s %12s\n','model','w0','t_rev (s)','w3(end)');
fprintf('%s\n',repmat('-',1,60));
for j = 1:2
    ax = nexttile(tl); hold(ax,'on');
    hh = gobjects(numel(w0List),1); nm = cell(numel(w0List),1);
    for k = 1:numel(w0List)
        P = P0; P.law = laws{j}; P.rho = 0;
        [T,Y,~] = run_case(t,P,-w0List(k),CFG);
        R = reversal(T,Y(:,3),-w0List(k));
        fprintf('%24s %8.1f %12s %12.3f\n', ...
                names{j}, w0List(k), fmt_t(R.t_cross), Y(end,3));
        s = 0.35 + 0.65*(k-1)/max(1,numel(w0List)-1);
        cc = cols{j}*s + (1-s)*[0.78 0.83 0.88];
        hh(k) = plot(ax,T,Y(:,3),'-','Color',cc,'LineWidth',1.5);
        nm{k} = sprintf('\\omega_0 = %.1f rad/s',w0List(k));
    end
    plot(ax,[t(1) t(end)],[0 0],':','Color',CFG.col.grey, ...
         'LineWidth',0.8,'HandleVisibility','off');
    hold(ax,'off'); style(ax,CFG);
    legend(ax,hh,nm,'Location','southeast');
    xlabel(ax,'$t$ (s)','Interpreter','latex');
    ylabel(ax,'$\omega_3$ (rad/s)','Interpreter','latex');
    title(ax,names{j});
end
savefig_(f,'compare_initial_spin',CFG);
end

% =========================================================================
% STUDY 6 - SOLVED NORMAL REACTION AGAINST N = m*g
% Holding N at m*g sets the vertical acceleration of G to zero while the
% height h(n) keeps oscillating, so the contact condition is violated. The
% error is small when friction is weak (little rocking) and grows towards
% the rolling limit, where the solved N converges to rolling and N = m*g
% does not.
% =========================================================================
function study_fixN(P0,CFG)
banner('FIXED N - solved normal reaction against N = m*g');
muList = [0.15 0.2 0.5 1.0 2.0];
t = linspace(0,CFG.tEnd,CFG.nSweep);

Pr = P0; Pr.law = 'rolling';
[Tr,Yr,~] = run_case(t,Pr,-CFG.w0,CFG);
Rr = reversal(Tr,Yr(:,3),-CFG.w0);

trS = nan(size(muList)); trF = nan(size(muList));
fprintf('%8s %16s %16s %10s\n','mu','t_rev solved N','t_rev N = mg','error');
fprintf('%s\n',repmat('-',1,54));
for k = 1:numel(muList)
    P = P0; P.law = 'coulomb'; P.rho = 0; P.mu = muList(k);
    P.fixN = false; [T,Y,~]  = run_case(t,P,-CFG.w0,CFG);
    P.fixN = true;  [~,Yf,~] = run_case(t,P,-CFG.w0,CFG);
    RS = reversal(T,Y(:,3),-CFG.w0);  trS(k) = RS.t_cross;
    RF = reversal(T,Yf(:,3),-CFG.w0); trF(k) = RF.t_cross;
    fprintf('%8.2f %16s %16s %9.1f%%\n',muList(k),fmt_t(trS(k)), ...
            fmt_t(trF(k)),100*(trF(k)-trS(k))/trS(k));
    if muList(k) == 1.0, YS = Y; YF = Yf; TS = T; end
end
fprintf('rolling benchmark: t_rev = %s s\n',fmt_t(Rr.t_cross));

f  = newfig('fixed N',[680 520],CFG);
tl = tiledlayout(f,2,1,'Padding','compact','TileSpacing','compact');

ax = nexttile(tl); hold(ax,'on');
h1 = plot(ax,Tr,Yr(:,3),'-','Color',CFG.col.ink,'LineWidth',2.2);
h2 = plot(ax,TS,YS(:,3),'-','Color',CFG.col.blue,'LineWidth',1.4);
h3 = plot(ax,TS,YF(:,3),'--','Color',CFG.col.verm,'LineWidth',1.4);
plot(ax,[t(1) t(end)],[0 0],':','Color',CFG.col.grey, ...
     'LineWidth',0.8,'HandleVisibility','off');
hold(ax,'off'); style(ax,CFG);
legend(ax,[h1 h2 h3],{'Rolling','Coulomb \mu = 1, solved N', ...
       'Coulomb \mu = 1, N = mg'},'Location','southeast');
xlabel(ax,'$t$ (s)','Interpreter','latex');
ylabel(ax,'$\omega_3$ (rad/s)','Interpreter','latex');
title(ax,'Spin reversal with solved and fixed normal reaction');

ax = nexttile(tl); hold(ax,'on');
hA = plot(ax,muList,trS,'o-','Color',CFG.col.blue,'LineWidth',1.6, ...
          'MarkerFaceColor','w','MarkerSize',6);
hB = plot(ax,muList,trF,'s--','Color',CFG.col.verm,'LineWidth',1.6, ...
          'MarkerFaceColor','w','MarkerSize',6);
hC = plot(ax,[muList(1) muList(end)],Rr.t_cross*[1 1],':', ...
          'Color',CFG.col.ink,'LineWidth',1.2);
hold(ax,'off'); style(ax,CFG);
legend(ax,[hA hB hC],{'solved N','N = mg','rolling'},'Location','northeast');
xlabel(ax,'$\mu$','Interpreter','latex');
ylabel(ax,'$t_{\mathrm{rev}}$ (s)','Interpreter','latex');
title(ax,'Reversal time against friction coefficient');
savefig_(f,'compare_fixedN',CFG);
end

% =========================================================================
% RUNNER
% Returns a 9-column Y for every model so that plotting, diagnostics and
% animation are identical whether the contact rolls or slips. For the
% rolling model the CM velocity is reconstructed from v_G = -w x r.
% =========================================================================
function [T,Y,D] = run_case(t,P,w3,CFG)
opts = odeset('AbsTol',CFG.absTol,'RelTol',CFG.relTol);
if strcmpi(P.law,'rolling')
    o6 = odeset('AbsTol',1e-12,'RelTol',1e-10);     % not stiff, ode45 fine
    [T,Y6] = ode45(@(tt,y) equns123_rolling(tt,y,P), t, ...
                   [0;0;w3;1e-6;1e-6;1e-6], o6);
    Y = [Y6, zeros(size(Y6,1),3)];
    for k = 1:size(Y6,1)
        w = Y6(k,1:3).'; ang = Y6(k,4:6).';
        [~,~,r,~,~] = kinematics(ang,w,P);
        Y(k,7:9) = (-cross(w,r)).';
    end
    D = diagnostics_rolling(T,Y6,P);
else
    % ode15s: the regularised Coulomb law has Jacobian entries of order
    % mu*N/eps (~5e4 N per m/s at eps = 3e-5), which is stiff. An explicit
    % solver such as ode45 will not finish in reasonable time.
    [T,Y] = ode15s(@(tt,y) equns123(tt,y,P), t, initcond(w3,P), opts);
    D = diagnostics(T,Y,P);
end
end

function y0 = initcond(w3,P)
% Start from a rolling state, v_P = v_G + w x r = 0, so that the large-mu
% limit can be compared like-for-like against the rolling benchmark.
ang = [1e-6; 1e-6; 1e-6];
w   = [0; 0; w3];
[~,~,r,~,~] = kinematics(ang,w,P);
y0  = [w; ang; -cross(w,r)];
end

% =========================================================================
% RIGHT-HAND SIDE - 9 first order ODEs
%   y = [w1 w2 w3 phi theta psi v1 v2 v3]'   (v in BODY components)
% =========================================================================
function dy = equns123(~,y,P)
w = y(1:3); ang = y(4:6); v = y(7:9);
[wd,vd] = contact_solve(w,ang,v,P);
theta = ang(2); psi = ang(3);
Phim  = [ cos(theta)*cos(psi)  sin(psi)  0;
         -cos(theta)*sin(psi)  cos(psi)  0;
          sin(theta)           0         1];
dy = [wd; Phim\w; vd];
end

% =========================================================================
% GEOMETRY AND ITS RATES
% =========================================================================
function [n,nd,r,rd,h] = kinematics(ang,w,P)
phi = ang(1); theta = ang(2); psi = ang(3);
a2v = P.a.^2;
n = [ sin(phi)*sin(psi) - cos(psi)*sin(theta)*cos(phi);
      cos(psi)*sin(phi) + sin(psi)*sin(theta)*cos(phi);
      cos(theta)*cos(phi) ];
nd = -cross(w,n);                        % since E3dot = 0 in space
S  = a2v.'*(n.^2);
sq = sqrt(S);
r  = -a2v.*n/sq;  r(3) = r(3) + P.b;     % r = r_{P/G}
Sd = a2v.'*(n.*nd);
rd = -a2v.*nd/sq + a2v.*n*Sd/S^1.5;      % b constant -> no term in rd(3)
h  = sq - P.b*n(3);                      % height of G above plane = -r.n
end

% =========================================================================
% CONTACT SOLVE - normal reaction from the contact condition
% =========================================================================
function [wd,vd,N,Ffr,Tsp,vp,gv] = contact_solve(w,ang,v,P)
[n,nd,r,rd,~] = kinematics(ang,w,P);

vp = v + cross(w,r);                     % slip velocity of P
vp = vp - (vp.'*n)*n;                    % tangential part only

d = cross(n,r);                          % zGddot = c + d.wdot
c = -(rd.'*nd) + r.'*cross(w,nd);

A0 = P.I\(-cross(w,P.I*w));              % gyroscopic contribution to wdot
B  = P.I\ cross(r,n);                    % contribution per unit N

sgn = 0; uhat = zeros(3,1); Fv = [];
switch lower(P.law)
    case 'coulomb'
        uhat = vp/sqrt(vp.'*vp + P.eps^2);            % regularised sign
        tau  = cross(r,-P.mu*uhat);                   % moment per unit N
        if P.rho > 0                                  % drilling friction
            sc   = P.rho*(w.'*n);
            sgn  = sc/sqrt(sc^2 + vp.'*vp + P.eps^2); % Contensou coupling
            tau  = tau - P.mu*P.rho*sgn*n;
        end
        A1 = P.I\tau;
    case 'viscous'
        Fv = -P.cv*vp;                                % independent of N ...
        A0 = P.I\(cross(r,Fv) - cross(w,P.I*w));      % ... so it goes in A0
        A1 = zeros(3,1);
    otherwise
        error('Unknown friction law: %s',P.law);
end

den = 1 - P.m*(d.'*B) - P.m*(d.'*A1);
N   = P.m*(P.g + c + d.'*A0)/den;
if P.fixN, N = P.m*P.g; end              % comparison only, see study_fixN

switch lower(P.law)
    case 'coulomb'
        Ffr = -P.mu*N*uhat;
        Tsp = zeros(3,1);
        if P.rho > 0, Tsp = -P.mu*N*P.rho*sgn*n; end
    case 'viscous'
        Ffr = Fv;
        Tsp = zeros(3,1);
end

wd = A0 + N*(B + A1);
fP = N*n + Ffr;
vd = fP/P.m - P.g*n - cross(w,v);        % body-frame derivative of v_G
gv = v.'*n + r.'*nd;                     % contact residual, exactly zero
vd = vd - P.kb*gv*n;                     % optional drift stabilisation
end

% =========================================================================
% ROLLING MODEL - six first-order ODEs with v_P = 0
% =========================================================================
function dy = equns123_rolling(~,y,P)
w = y(1:3); ang = y(4:6);
[n,~,r,rd,~] = kinematics(ang,w,P);
Xm = [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0];
Wm = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
Ir    = P.I - P.m*(Xm*Xm);
rhsvc = -Wm*P.I*w - P.m*(Xm*(Wm*(Wm*r))) - P.m*(Xm*(Wm*rd)) + P.m*P.g*(Xm*n);
wd    = Ir\rhsvc;
theta = ang(2); psi = ang(3);
Phim  = [ cos(theta)*cos(psi)  sin(psi)  0;
         -cos(theta)*sin(psi)  cos(psi)  0;
          sin(theta)           0         1];
dy = [wd; Phim\w];
end

% =========================================================================
% DIAGNOSTICS - contact force, slip, energy and residual at every output point
% =========================================================================
function D = diagnostics(T,Y,P)
step  = 1;          % every output point: the peaks of N last only a few ms
D.idx = 1:step:size(Y,1);
D.T   = T(D.idx);
nT    = numel(D.idx);
[D.N,D.F,D.Tsp,D.slip,D.E,D.gv] = deal(zeros(nT,1));
for kk = 1:nT
    k = D.idx(kk);
    w = Y(k,1:3).'; ang = Y(k,4:6).'; v = Y(k,7:9).';
    [~,~,N,Ffr,Tsp,vp,gv] = contact_solve(w,ang,v,P);
    [~,~,~,~,h] = kinematics(ang,w,P);
    D.N(kk) = N;          D.F(kk)  = norm(Ffr);  D.Tsp(kk) = norm(Tsp);
    D.slip(kk) = norm(vp); D.gv(kk) = gv;
    D.E(kk) = 0.5*w.'*P.I*w + 0.5*P.m*(v.'*v) + P.m*P.g*h;
end
end

function D = diagnostics_rolling(T,Y6,P)
% The contact force is recovered a posteriori from f_P = m*a_G + m*g*E3,
% the quantity eliminated in the rolling formulation. This makes the
% normal reaction of the rolling solution observable.
step  = 1;
D.idx = 1:step:size(Y6,1);
D.T   = T(D.idx);
nT    = numel(D.idx);
[D.N,D.F,D.Tsp,D.slip,D.E,D.gv] = deal(zeros(nT,1));
for kk = 1:nT
    k = D.idx(kk);
    w = Y6(k,1:3).'; ang = Y6(k,4:6).';
    [n,~,r,rd,h] = kinematics(ang,w,P);
    dy = equns123_rolling(0,[w;ang],P);
    wd = dy(1:3);
    aG = -cross(wd,r) - cross(w,cross(w,r)) - cross(w,rd);
    fP = P.m*aG + P.m*P.g*n;
    v  = -cross(w,r);
    D.N(kk) = fP.'*n;
    D.F(kk) = norm(fP - (fP.'*n)*n);
    D.E(kk) = 0.5*w.'*P.I*w + 0.5*P.m*(v.'*v) + P.m*P.g*h;
end
end

% =========================================================================
% REVERSAL DETECTION
% A zero crossing alone is NOT a reversal. Drilling friction can spin the
% body down through zero without reversing it at all. A true reversal also
% requires the spin to settle at a finite magnitude in the opposite sense.
% =========================================================================
function R = reversal(T,w3,w30)
opp = -sign(w30);
idx = find(sign(w3) == opp & T > 0.05, 1, 'first');
if isempty(idx)
    R.t_cross = NaN; R.settled = 0; R.true = false;
else
    R.t_cross = T(idx);
    R.settled = max(abs(w3(idx:end)));
    R.true    = R.settled > 0.05;   % rad/s; Coulomb friction recovers ~0.3 rad/s,
                                    % so a higher threshold would misclassify it
end
end

function s = fmt_t(t)
if isnan(t)
    s = 'none';
else
    s = sprintf('%.2f',t);
end
end

% =========================================================================
% CONSOLE TABLE
% =========================================================================
function hdr_table()
fprintf('\n%-30s %-14s %-9s %-9s %-18s %-11s %s\n', ...
        'model','sense','t_rev','w3(end)','N range (N)','dE (J)','checks');
fprintf('%s\n',repmat('-',1,112));
end

function row_table(name,sense,T,Y,D,w30)
R = reversal(T,Y(:,3),w30);
if isnan(R.t_cross)
    tr = 'none';
elseif ~R.true
    tr = sprintf('%.2f*',R.t_cross);      % crossing, but not a reversal
else
    tr = sprintf('%.2f',R.t_cross);
end
if all(diff(D.E) < 1e-9), eok = 'E-ok'; else, eok = 'E-BAD'; end
if min(D.N) > 0,          nok = 'N-ok'; else, nok = 'N-BAD'; end
fprintf('%-30s %-14s %-9s %-9.3f [%7.2f,%8.2f] %-11.2e gv=%.0e %s %s\n', ...
        name, sense, tr, Y(end,3), min(D.N), max(D.N), ...
        D.E(end)-D.E(1), max(abs(D.gv)), eok, nok);
if ~isnan(R.t_cross) && ~R.true
    fprintf('%30s %s\n','','* zero crossing is SPIN-DOWN, not a reversal');
end
end

function banner(txt)
fprintf('\n%s\n%s\n%s\n',repmat('=',1,72),txt,repmat('=',1,72));
end

% =========================================================================
% PLOTTING
% -------------------------------------------------------------------------
%  * Legends are built from explicit graphics handles, legend(ax,h,names),
%    so that reference lines added to an axes cannot shift the mapping
%    between legend entries and curves.
%  * omega_1 and omega_2 oscillate at about 11 Hz, so a raw trace over 8 s
%    is a solid band. The moving-window envelope is drawn over a faded raw
%    trace so that the modulation is visible.
%  * Okabe-Ito palette (legible in greyscale), tiled layouts, 300 dpi export.
% =========================================================================

function f = newfig(name,sz,CFG)
f = figure('Name',name,'Color','w','Units','pixels', ...
           'Position',[80 80 sz(1) sz(2)]);
set(f,'DefaultAxesFontName',CFG.font,'DefaultAxesFontSize',CFG.fs, ...
      'DefaultTextFontName',CFG.font);
end

function style(ax,CFG)
set(ax,'FontName',CFG.font,'FontSize',CFG.fs,'LineWidth',0.8, ...
       'TickDir','out','Box','off','Layer','top', ...
       'XColor',[0.25 0.28 0.32],'YColor',[0.25 0.28 0.32]);
grid(ax,'on');
set(ax,'GridColor',[0.60 0.65 0.70],'GridAlpha',0.25);
% Write tick labels in full. Otherwise MATLAB may move a common factor such
% as x10^-3 to the top of the axis, where the panel title hides it and the
% axis appears to read 1000 times too large.
try
    ax.YAxis.Exponent = 0;
    ytickformat(ax,'%g');
catch
end
end

function savefig_(f,tag,CFG)
if ~CFG.export, return; end
if ~isfolder(CFG.figdir), mkdir(CFG.figdir); end
fn = fullfile(CFG.figdir,[tag '.png']);
try
    exportgraphics(f,fn,'Resolution',CFG.dpi,'BackgroundColor','white');
catch ME
    warning('exportgraphics failed (%s); using print instead.',ME.message);
    print(f,fn,'-dpng',sprintf('-r%d',CFG.dpi));
end
end

function d = prepare_figdir(name)
% Resolve the output folder to an ABSOLUTE path next to this .m file, so the
% result does not depend on MATLAB's current folder. (exist(name,'dir') is
% not used: it also returns true for a folder of that name anywhere on the
% MATLAB path, in which case nothing would be created in the current folder.)
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
d = fullfile(here,name);
if ~isfolder(d)
    [ok,msg] = mkdir(d);
    if ~ok                                   % e.g. a read-only folder
        warning('Cannot create %s (%s).',d,msg);
        d = fullfile(tempdir,'rattleback_figures');
        if ~isfolder(d), mkdir(d); end
    end
end
fprintf('Figures will be written to:\n  %s\n',d);
end

function [lo,hi] = envel(y,win)
% envelope of a fast oscillation, so the figure shows the modulation
% rather than a solid band
win = max(3,round(win));
hi  = movmax(y,win);
lo  = movmin(y,win);
end

% ---- per-model angular velocity figure ----------------------------------
function fig_omega(T,Y1,Y2,M,CFG)
f   = newfig(['omega: ' M.tag],[720 500],CFG);
tl  = tiledlayout(f,2,1,'Padding','compact','TileSpacing','compact');
win = round(0.15/(T(2)-T(1)));     % 1.5 periods of the 11 Hz rocking
cols = {CFG.col.blue, CFG.col.verm, CFG.col.green};
labs = {'\omega_1','\omega_2','\omega_3'};
Ys   = {Y1,Y2};
ttl  = {'preferred sense','non-preferred sense'};

for s = 1:2
    ax = nexttile(tl); hold(ax,'on');
    Y = Ys{s}; hh = gobjects(3,1);
    for j = 1:3
        c = cols{j};
        if j < 3
            [lo,hi] = envel(Y(:,j),win);
            fill(ax,[T; flipud(T)],[hi; flipud(lo)],c, ...
                 'FaceAlpha',0.16,'EdgeColor','none','HandleVisibility','off');
            plot(ax,T,Y(:,j),'-','Color',[c 0.30],'LineWidth',0.5, ...
                 'HandleVisibility','off');
            hh(j) = plot(ax,T,hi,'-','Color',c,'LineWidth',1.1);
            plot(ax,T,lo,'-','Color',c,'LineWidth',1.1,'HandleVisibility','off');
        else
            hh(j) = plot(ax,T,Y(:,j),'-','Color',c,'LineWidth',1.8);
        end
    end
    plot(ax,[T(1) T(end)],[0 0],':','Color',CFG.col.grey, ...
         'LineWidth',0.8,'HandleVisibility','off');
    hold(ax,'off'); style(ax,CFG);
    legend(ax,hh,labs,'Location','best','Orientation','horizontal');
    xlabel(ax,'$t$ (s)','Interpreter','latex');
    ylabel(ax,'$\omega$ (rad/s)','Interpreter','latex');
    title(ax,[M.name ' - ' ttl{s}]);
end
savefig_(f,['omega_' M.tag],CFG);
end

% ---- per-model nutation figure ------------------------------------------
function fig_theta(T,Y1,Y2,M,CFG)
f   = newfig(['theta: ' M.tag],[720 450],CFG);
tl  = tiledlayout(f,2,1,'Padding','compact','TileSpacing','compact');
win = round(0.15/(T(2)-T(1)));     % 1.5 periods of the 11 Hz rocking
Ys  = {Y1,Y2};
ttl = {'preferred sense','non-preferred sense'};
cs  = {CFG.col.pref, CFG.col.nonpref};

for s = 1:2
    ax = nexttile(tl); hold(ax,'on');
    th = Ys{s}(:,5)*180/pi; c = cs{s};
    [lo,hi] = envel(th,win);
    fill(ax,[T; flipud(T)],[hi; flipud(lo)],c, ...
         'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
    plot(ax,T,th,'-','Color',[c 0.40],'LineWidth',0.5,'HandleVisibility','off');
    h1 = plot(ax,T,hi,'-','Color',c,'LineWidth',1.2);
    plot(ax,T,lo,'-','Color',c,'LineWidth',1.2,'HandleVisibility','off');
    plot(ax,[T(1) T(end)],[0 0],':','Color',CFG.col.grey, ...
         'LineWidth',0.8,'HandleVisibility','off');
    hold(ax,'off'); style(ax,CFG);
    legend(ax,h1,{'envelope of \theta'},'Location','best');
    xlabel(ax,'$t$ (s)','Interpreter','latex');
    ylabel(ax,'$\theta$ (deg)','Interpreter','latex');
    title(ax,[M.name ' - ' ttl{s}]);
end
savefig_(f,['theta_' M.tag],CFG);
end

% ---- per-model contact diagnostics --------------------------------------
function fig_contact(D1,D2,M,CFG)
mg = 9.81;
f  = newfig(['contact: ' M.tag],[800 560],CFG);
tl = tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');

Ap = {D1.N/mg, D1.F, D1.slip, D1.E};
An = {D2.N/mg, D2.F, D2.slip, D2.E};
yl = {'$N/mg$','$|\mathbf{F}|$ (N)','$|\mathbf{v}_P|$ (m/s)','$E$ (J)'};
tt = {'normal force','tangential force at P', ...
      'slip velocity at P','total energy (must not increase)'};

for k = 1:4
    ax = nexttile(tl); hold(ax,'on');
    h1 = plot(ax,D1.T,Ap{k},'-','Color',CFG.col.pref,'LineWidth',1.0);
    h2 = plot(ax,D2.T,An{k},'-','Color',CFG.col.nonpref,'LineWidth',1.0);
    if k == 1
        h3 = plot(ax,[D1.T(1) D1.T(end)],[1 1],'--', ...
                  'Color',CFG.col.ink,'LineWidth',1.0);
        legend(ax,[h1 h2 h3],{'preferred','non-preferred','N = mg'}, ...
               'Location','best');
    else
        legend(ax,[h1 h2],{'preferred','non-preferred'},'Location','best');
    end
    hold(ax,'off'); style(ax,CFG);
    if k == 4
        % the rolling energy is constant to ~1e-11 J; give the axis a
        % sensible minimum span so round-off is not magnified into noise
        Eall = [D1.E; D2.E]; lo = min(Eall); hi = max(Eall);
        span = max(hi - lo, 0.15*hi); mid = 0.5*(lo + hi);
        if hi - lo < 0.15*hi, ylim(ax, mid + 0.5*span*[-1 1]);
        else, ylim(ax, [lo - 0.05*span, hi + 0.05*span]); end
    end
    xlabel(ax,'$t$ (s)','Interpreter','latex');
    ylabel(ax,yl{k},'Interpreter','latex');
    title(ax,tt{k});
end
title(tl,M.name,'FontWeight','normal','FontName',CFG.font);
savefig_(f,['contact_' M.tag],CFG);
end

% ---- cross-model contact comparison -------------------------------------
function fig_compare_contact(MDL,CFG)
mg = 9.81; nM = numel(MDL);
f  = newfig('friction model comparison',[900 620],CFG);
tl = tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
nm = cellfun(@(s) s.name, MDL,'UniformOutput',false);

% spin reversal
ax = nexttile(tl); hold(ax,'on'); hh = gobjects(nM,1);
for k = 1:nM
    hh(k) = plot(ax,MDL{k}.T,MDL{k}.Y2(:,3),'-', ...
                 'Color',MDL{k}.col,'LineWidth',1.6);
end
plot(ax,[MDL{1}.T(1) MDL{1}.T(end)],[0 0],':','Color',CFG.col.grey, ...
     'LineWidth',0.8,'HandleVisibility','off');
hold(ax,'off'); style(ax,CFG);
legend(ax,hh,nm,'Location','southeast');
xlabel(ax,'$t$ (s)','Interpreter','latex');
ylabel(ax,'$\omega_3$ (rad/s)','Interpreter','latex');
title(ax,'spin reversal');

% tangential force
ax = nexttile(tl); hold(ax,'on'); hh = gobjects(nM,1);
for k = 1:nM
    hh(k) = plot(ax,MDL{k}.D2.T,MDL{k}.D2.F,'-', ...
                 'Color',MDL{k}.col,'LineWidth',1.0);
end
hold(ax,'off'); style(ax,CFG);
legend(ax,hh,nm,'Location','northeast');
xlabel(ax,'$t$ (s)','Interpreter','latex');
ylabel(ax,'$|\mathbf{F}|$ (N)','Interpreter','latex');
title(ax,'tangential force at P');

% normal force
ax = nexttile(tl); hold(ax,'on'); hh = gobjects(nM+1,1);
for k = 1:nM
    hh(k) = plot(ax,MDL{k}.D2.T,MDL{k}.D2.N/mg,'-', ...
                 'Color',MDL{k}.col,'LineWidth',1.0);
end
hh(nM+1) = plot(ax,[MDL{1}.T(1) MDL{1}.T(end)],[1 1],'--', ...
                'Color',[0.35 0.38 0.42],'LineWidth',1.0);
hold(ax,'off'); style(ax,CFG);
legend(ax,hh,[nm, {'N = mg'}],'Location','northeast');
xlabel(ax,'$t$ (s)','Interpreter','latex');
ylabel(ax,'$N/mg$','Interpreter','latex');
title(ax,'normal force (N = mg would be the dashed line)');

% slip velocity
ax = nexttile(tl); hold(ax,'on'); hh = gobjects(nM-1,1);
for k = 2:nM
    hh(k-1) = plot(ax,MDL{k}.D2.T,MDL{k}.D2.slip,'-', ...
                   'Color',MDL{k}.col,'LineWidth',1.0);
end
hold(ax,'off'); style(ax,CFG);
legend(ax,hh,nm(2:end),'Location','northeast');
xlabel(ax,'$t$ (s)','Interpreter','latex');
ylabel(ax,'$|\mathbf{v}_P|$ (m/s)','Interpreter','latex');
title(ax,'slip velocity (identically zero for rolling)');

title(tl,'Non-preferred sense','FontWeight','normal','FontName',CFG.font);
savefig_(f,'comparison_contact',CFG);
end

% ---- cross-model angular velocity comparison ----------------------------
function fig_compare_omega(MDL,CFG)
nM = numel(MDL);
f  = newfig('angular velocity comparison',[900 680],CFG);
tl = tiledlayout(f,3,1,'Padding','compact','TileSpacing','compact');
nm  = cellfun(@(s) s.name, MDL,'UniformOutput',false);
lab = {'$\omega_1$ (rad/s)','$\omega_2$ (rad/s)','$\omega_3$ (rad/s)'};
win = round(0.15/(MDL{1}.T(2)-MDL{1}.T(1)));

for j = 1:3
    ax = nexttile(tl); hold(ax,'on'); hh = gobjects(nM,1);
    for k = 1:nM
        y = MDL{k}.Y2(:,j); T = MDL{k}.T; c = MDL{k}.col;
        if j < 3
            [~,hi] = envel(y,win);            % upper envelope only:
            hh(k) = plot(ax,T,hi,'-','Color',c,'LineWidth',1.3);
            plot(ax,T,-hi,'-','Color',[c 0.30],'LineWidth',0.7, ...
                 'HandleVisibility','off');
        else
            hh(k) = plot(ax,T,y,'-','Color',c,'LineWidth',1.5);
        end
    end
    plot(ax,[MDL{1}.T(1) MDL{1}.T(end)],[0 0],':','Color',CFG.col.grey, ...
         'LineWidth',0.8,'HandleVisibility','off');
    hold(ax,'off'); style(ax,CFG);
    if j == 1
        legend(ax,hh,nm,'Location','northeast');
        title(ax,'envelope of the rattling (non-preferred sense)');
    end
    xlabel(ax,'$t$ (s)','Interpreter','latex');
    ylabel(ax,lab{j},'Interpreter','latex');
end
savefig_(f,'comparison_omega',CFG);
end

% =========================================================================
% ANIMATION
%   * rotation assembled as Q = R1(phi)*R2(theta)*R3(psi), the 1-2-3
%     sequence. Q'*E3 reproduces n of kinematics() to machine precision
%   * G placed at its true height h(n), and the red marker is the COMPUTED
%     contact point, so if it stays on the table the constraint is holding
%   * CM ground track integrated from the solution, so slipping shows as
%     visible translation and rolling traces a rolling path
%   * true semi-ellipsoid at the actual semi-axes, lower half plus flat cap
% =========================================================================
function z = animate(T,Y,P,tag,CFG)
z = 0; nframes = CFG.animFrames;

nT = size(Y,1); vx = zeros(nT,1); vy = zeros(nT,1);
for k = 1:nT
    vs = rotmat(Y(k,4:6))*Y(k,7:9).';
    vx(k) = vs(1); vy(k) = vs(2);
end
xg = cumtrapz(T,vx); yg = cumtrapz(T,vy);

nu = 72; nv = 28;
u = linspace(0,2*pi,nu); v = linspace(-pi/2,0,nv);
[U,V] = meshgrid(u,v);
Xb = P.a1*cos(V).*cos(U);
Yb = P.a2*cos(V).*sin(U);
Zb = P.a3*sin(V) + P.b;                % x3 = xbreve3 + b
sr = linspace(0,1,8);
[S,Uc] = meshgrid(sr,u);
Xc = P.a1*S.*cos(Uc); Yc = P.a2*S.*sin(Uc); Zc = P.b*ones(size(S));

idx = round(linspace(1,nT,nframes));
L   = 1.7*P.a1;
fig = figure('Color','w','Name',['animation: ' tag], ...
             'Position',[80 80 720 540]);
map = []; mov = [];

for jj = 1:numel(idx)
    k   = idx(jj);
    ang = Y(k,4:6).'; w = Y(k,1:3).';
    Q   = rotmat(ang);
    [~,~,r,~,h] = kinematics(ang,w,P);
    o = [xg(k); yg(k); h];
    [Xs,Ys,Zs] = xform(Q,o,Xb,Yb,Zb);
    [Xk,Yk,Zk] = xform(Q,o,Xc,Yc,Zc);
    rp = Q*r + o;

    clf(fig);
    ax = axes('Parent',fig); hold(ax,'on');
    patch(ax,[-1 1 1 -1]*L*3,[-1 -1 1 1]*L*3,[0 0 0 0], ...
          [0.91 0.94 0.96],'EdgeColor','none');
    gx = round(xg(k)/0.05)*0.05; gy = round(yg(k)/0.05)*0.05;
    for gl = -6:6
        plot3(ax,[gx gx]+gl*0.05, gy+[-6 6]*0.05, [0 0], ...
              '-','Color',[0.80 0.85 0.89],'LineWidth',0.5);
        plot3(ax,gx+[-6 6]*0.05, [gy gy]+gl*0.05, [0 0], ...
              '-','Color',[0.80 0.85 0.89],'LineWidth',0.5);
    end
    surf(ax,Xs,Ys,Zs,'EdgeColor','none','FaceColor',[0.70 0.44 0.20], ...
         'FaceLighting','gouraud','AmbientStrength',0.45, ...
         'SpecularStrength',0.30);
    surf(ax,Xk,Yk,Zk,'EdgeColor','none','FaceColor',[0.85 0.62 0.35], ...
         'FaceLighting','gouraud','AmbientStrength',0.50);
    plot3(ax,xg(1:k),yg(1:k),zeros(k,1),'-', ...
          'Color',[0.35 0.42 0.50],'LineWidth',1.2);
    plot3(ax,rp(1),rp(2),rp(3),'o','MarkerFaceColor',[0.83 0.25 0.29], ...
          'MarkerEdgeColor','none','MarkerSize',6);
    camlight headlight; material dull;
    hold(ax,'off');

    daspect(ax,[1 1 1]);
    xlim(ax,xg(k)+[-L L]); ylim(ax,yg(k)+[-L L]); zlim(ax,[0 1.3*L]);
    view(ax,38,20); axis(ax,'off');
    title(ax,sprintf('%s      t = %5.2f s      \\omega_3 = %+6.2f rad/s', ...
          strrep(tag,'_','\_'), T(k), Y(k,3)), ...
          'FontName',CFG.font,'FontSize',11,'FontWeight','normal');

    fr = getframe(fig);
    if jj == 1
        [mov(:,:,1,1), map] = rgb2ind(fr.cdata,256,'nodither');
    else
        mov(:,:,1,jj) = rgb2ind(fr.cdata,map,'nodither');
    end
end
for i = 1:12
    mov(:,:,1,numel(idx)+i) = mov(:,:,1,numel(idx));
end

fn = ['animation_' tag '.gif'];
imwrite(mov,map,fn,'DelayTime',0.04,'LoopCount',inf);
fprintf('wrote %s\n',fn);
end

% ---- body-to-space rotation, 1-2-3 sequence --------------------------------
function Q = rotmat(ang)
phi = ang(1); th = ang(2); ps = ang(3);
R1 = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
R2 = [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
R3 = [cos(ps) -sin(ps) 0; sin(ps) cos(ps) 0; 0 0 1];
Q  = R1*R2*R3;
end

function [Xs,Ys,Zs] = xform(Q,o,Xb,Yb,Zb)
Pb = [Xb(:).'; Yb(:).'; Zb(:).'];
Pw = Q*Pb + o;
Xs = reshape(Pw(1,:),size(Xb));
Ys = reshape(Pw(2,:),size(Xb));
Zs = reshape(Pw(3,:),size(Xb));
end
