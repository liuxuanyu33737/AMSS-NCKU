# MPI 进程 829 的 gprof 性能报告

- 数据文件：`gmon.out.829`
- 函数数量：131
- 分析工具：GNU gprof

|排名|占比(%)|累计占比(%)|自身时间(s)|调用次数|函数名称|
|---:|---:|---:|---:|---:|---|
|1|57.26|57.26|222.77|5457|`compute_rhs_bssn_`|
|2|9.02|66.28|35.09|130968|`lopsided_`|
|3|7.26|73.54|28.23|60027|`fdderivs_`|
|4|5.99|79.53|23.31|130968|`kodis_`|
|5|5.31|84.84|20.65|91237|`fderivs_`|
|6|5.11|89.95|19.88|203393|`prolong3_`|
|7|2.27|92.22|8.83|5280|`enforce_ga_`|
|8|2.25|94.47|8.77|126720|`rungekutta4_rout_`|
|9|1.78|96.25|6.91|744441|`symmetry_bd_`|
|10|0.86|97.11|3.35|125928|`restrict3_`|
|11|0.79|97.90|3.08||`_init`|
|12|0.53|98.43|2.06|15936|`average2_`|
|13|0.29|98.72|1.14|521994240|`misc::fact(int)`|
|14|0.26|98.98|1.00|240|`admmass_bssn_`|
|15|0.20|99.18|0.77|20|`getnp4_`|
|16|0.19|99.37|0.74|1504474|`copy_`|
|17|0.11|99.48|0.41|240|`surface_integral::surf_Wave(double, int, cgh*, var*, var*, int, int, int, double*, double*, monitor*)`|
|18|0.08|99.56|0.32|1760|`Patch::Interp_Points(MyList<var>*, int, double**, double*, int)`|
|19|0.06|99.62|0.25|563166|`Ansorg::interpolate_tri_bar(double, double, double, int, int, int, double*, double*, double*, double*)`|
|20|0.06|99.68|0.22|46448640|`misc::Wigner_d_function(int, int, int, double)`|
|21|0.05|99.73|0.21|5280|`lowerboundset_`|
|22|0.05|99.78|0.18|2297060|`polint_`|
|23|0.03|99.81|0.12|1320|`bssn_class::Step(int, int)`|
|24|0.03|99.84|0.11|869188|`monitor::monitor(char const*, int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)`|
|25|0.03|99.87|0.11|124800|`sommerfeld_rout_`|
|26|0.03|99.90|0.10|29043|`Parallel::transfer(MyList<Parallel::gridseg>**, MyList<Parallel::gridseg>**, MyList<var>*, MyList<var>*, int)`|
|27|0.03|99.93|0.10|1920|`sommerfeld_routbam_`|
|28|0.02|99.95|0.06|1260|`l2normhelper_`|
|29|0.02|99.97|0.06|240|`surface_integral::surf_MassPAng(double, int, cgh*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, double*, monitor*)`|
|30|0.01|99.98|0.04|9|`get_ansorg_nbhs_`|
|31|0.01|99.99|0.02|1689498|`Ansorg::find_point_bisection(double, int, double*, int)`|
|32|0.01|100.00|0.02|231952|`Parallel::build_gstl(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>**, MyList<Parallel::gridseg>**)`|
|33|0.01|100.01|0.02|53420|`decide3d_`|
|34|0.01|100.02|0.02|53420|`polin3_`|
|35|0.01|100.03|0.02|19327|`Parallel::gsl_subtract(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*)`|
|36|0.01|100.04|0.02|1093|`Parallel::prepare_inter_time_level(Patch*, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, int)`|
|37|0.01|100.05|0.02|120|`average_`|
|38|0.00|100.05|0.01|13216975|`Block::getdX(int)`|
|39|0.00|100.05|0.01|150168|`Parallel::gs_subtract(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*)`|
|40|0.00|100.05|0.01|65696|`Parallel::build_bulk_gsl(Block*, Patch*)`|
|41|0.00|100.05|0.01|42240|`Block::swapList(MyList<var>*, MyList<var>*, int)`|
|42|0.00|100.05|0.01|21000|`cgh::Interp_One_Point(MyList<var>*, double*, double*, int)`|
|43|0.00|100.05|0.01|1406|`Parallel::Restrict(MyList<Patch>*, MyList<Patch>*, MyList<var>*, MyList<var>*, int)`|
|44|0.00|100.05|0.01|12|`Parallel::Dump_Data(Patch*, MyList<var>*, char*, double, double, int)`|
|45|0.00|100.05|0.00|701766|`Patch::getdX(int)`|
|46|0.00|100.05|0.00|563166|`Ansorg::xyz_to_ABp(double, double, double, double*, double*, double*)`|
|47|0.00|100.05|0.00|563166|`Ansorg::ps_u_at_xyz(double, double, double)`|
|48|0.00|100.05|0.00|230517|`Patch::Interp_ONE_Point(MyList<var>*, double*, double*, int)`|
|49|0.00|100.05|0.00|140352|`Parallel::build_owned_gsl2(Patch*, int)`|
|50|0.00|100.05|0.00|120144|`Parallel::build_owned_gsl0(Patch*, int)`|
|51|0.00|100.05|0.00|53420|`global_interp_`|
|52|0.00|100.05|0.00|35552|`Parallel::build_owned_gsl4(Patch*, int, int)`|
|53|0.00|100.05|0.00|21821|`Parallel::build_complete_gsl(Patch*)`|
|54|0.00|100.05|0.00|19278|`Parallel::build_buffer_gsl(Patch*)`|
|55|0.00|100.05|0.00|14923|`Parallel::build_ghost_gsl(Patch*)`|
|56|0.00|100.05|0.00|14923|`Parallel::Sync(Patch*, MyList<var>*, int)`|
|57|0.00|100.05|0.00|8212|`Parallel::Sync(MyList<Patch>*, MyList<var>*, int)`|
|58|0.00|100.05|0.00|4355|`Parallel::OutBdLow2Hi(Patch*, Patch*, MyList<var>*, MyList<var>*, int)`|
|59|0.00|100.05|0.00|1939|`misc::parse_parts(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, int&)`|
|60|0.00|100.05|0.00|1321|`setpbh(int, double**, double*, int)`|
|61|0.00|100.05|0.00|1308|`bssn_class::Constraint_Out()`|
|62|0.00|100.05|0.00|1280|`Parallel::PatList_Interp_Points(MyList<Patch>*, MyList<var>*, int, double**, double*, int)`|
|63|0.00|100.05|0.00|1260|`Parallel::L2Norm(Patch*, var*)`|
|64|0.00|100.05|0.00|700|`cgh::Regrid_Onelevel(int, int, int, double**, double**, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, bool, monitor*)`|
|65|0.00|100.05|0.00|640|`bssn_class::compute_Porg_rhs(double**, double**, var*, var*, var*, int)`|
|66|0.00|100.05|0.00|639|`misc::parse_parts(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, int&, int&, int&)`|
|67|0.00|100.05|0.00|464|`Block::Block(int, int*, double*, int, int, int, int, int)`|
|68|0.00|100.05|0.00|464|`Block::~Block()`|
|69|0.00|100.05|0.00|440|`monitor::writefile(double, int, double*)`|
|70|0.00|100.05|0.00|240|`monitor::writefile(double, int, double*, double*)`|
|71|0.00|100.05|0.00|167|`var::var(char const*, int, double, double, double)`|
|72|0.00|100.05|0.00|167|`var::~var()`|
|73|0.00|100.05|0.00|156|`Parallel::writefile(double, int, int, int, double, double, double, double, double, double, char*, double*)`|
|74|0.00|100.05|0.00|107|`Patch::Patch(int, int*, double*, int, bool, int)`|
|75|0.00|100.05|0.00|107|`Patch::~Patch()`|
|76|0.00|100.05|0.00|107|`Parallel::partition3(int*, int, int*, int, int*)`|
|77|0.00|100.05|0.00|65|`frame_dummy`|
|78|0.00|100.05|0.00|58|`cgh::construct_patchlist(int, int)`|
|79|0.00|100.05|0.00|58|`Parallel::KillBlocks(MyList<Patch>*)`|
|80|0.00|100.05|0.00|58|`Parallel::distribute(MyList<Patch>*, int, int, int, bool, int)`|
|81|0.00|100.05|0.00|58|`Parallel::add_ghost_touch(MyList<Parallel::gridseg>*&)`|
|82|0.00|100.05|0.00|58|`Parallel::cut_gsl(MyList<Parallel::gridseg>*&)`|
|83|0.00|100.05|0.00|58|`Parallel::merge_gsl(MyList<Parallel::gridseg>*&, double)`|
|84|0.00|100.05|0.00|53|`Parallel::checkgsl(MyList<Parallel::gridseg>*, bool)`|
|85|0.00|100.05|0.00|49|`cgh::recompose_cgh_Onelevel(int, int, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, int, bool)`|
|86|0.00|100.05|0.00|49|`Parallel::fill_level_data(MyList<Patch>*, MyList<Patch>*, MyList<Patch>*, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, int, bool, bool)`|
|87|0.00|100.05|0.00|49|`Parallel::build_complete_gsl_virtual(MyList<Patch>*)`|
|88|0.00|100.05|0.00|49|`Parallel::cut_gs(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*&)`|
|89|0.00|100.05|0.00|21|`bssn_class::Interp_Constraint(bool)`|
|90|0.00|100.05|0.00|21|`bssn_class::Compute_Constraint()`|
|91|0.00|100.05|0.00|21|`cgh::cgh(int, int, int, char*, int, monitor*)`|
|92|0.00|100.05|0.00|20|`bssn_class::Compute_Psi4(int)`|
|93|0.00|100.05|0.00|20|`bssn_class::RecursiveStep(int)`|
|94|0.00|100.05|0.00|20|`bssn_class::check_Stdin_Abort()`|
|95|0.00|100.05|0.00|20|`perf::MemoryUsage(unsigned long*, unsigned long*, unsigned long*, unsigned long*, unsigned long*, unsigned long*, int)`|
|96|0.00|100.05|0.00|20|`perf::sample_mem_usage(int)`|
|97|0.00|100.05|0.00|15|`var::setpropspeed(double)`|
|98|0.00|100.05|0.00|15|`std::_Rb_tree<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::_Select1st<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >::find(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&)`|
|99|0.00|100.05|0.00|13|`Parallel::aligncheck(double*, double*, int, double*, int*)`|
|100|0.00|100.05|0.00|12|`Patch::checkPatch(bool)`|
|101|0.00|100.05|0.00|9|`Parallel::Dump_Data(MyList<Patch>*, MyList<var>*, char*, double, double)`|
|102|0.00|100.05|0.00|7|`Parallel::prepare_inter_time_level(Patch*, MyList<var>*, MyList<var>*, MyList<var>*, int)`|
|103|0.00|100.05|0.00|5|`monitor::~monitor()`|
|104|0.00|100.05|0.00|2|`checkpoint::addvariablelist(MyList<var>*)`|
|105|0.00|100.05|0.00|2|`misc::inversearray(double*, int)`|
|106|0.00|100.05|0.00|2|`std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >::~pair()`|
|107|0.00|100.05|0.00|2|`std::pair<std::_Rb_tree_iterator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, bool> std::_Rb_tree<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::_Select1st<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >::_M_insert_unique<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >(std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >&&)`|
|108|0.00|100.05|0.00|1|`bssn_class::Initialize()`|
|109|0.00|100.05|0.00|1|`bssn_class::Read_Ansorg()`|
|110|0.00|100.05|0.00|1|`bssn_class::Setup_Black_Hole_position()`|
|111|0.00|100.05|0.00|1|`bssn_class::Evolve(int)`|
|112|0.00|100.05|0.00|1|`bssn_class::bssn_class(double, double, double, double, double, double, double, int, int, char*, double, double, double, int, int, int, double, double)`|
|113|0.00|100.05|0.00|1|`bssn_class::~bssn_class()`|
|114|0.00|100.05|0.00|1|`checkpoint::checkpoint(bool, char const*, int)`|
|115|0.00|100.05|0.00|1|`checkpoint::~checkpoint()`|
|116|0.00|100.05|0.00|1|`surface_integral::surface_integral(int)`|
|117|0.00|100.05|0.00|1|`surface_integral::~surface_integral()`|
|118|0.00|100.05|0.00|1|`cgh::compose_cgh(int)`|
|119|0.00|100.05|0.00|1|`cgh::settrfls(int)`|
|120|0.00|100.05|0.00|1|`cgh::read_bbox(int, char*)`|
|121|0.00|100.05|0.00|1|`cgh::sethandle(monitor*)`|
|122|0.00|100.05|0.00|1|`cgh::~cgh()`|
|123|0.00|100.05|0.00|1|`misc::gaulegf(double, double, double*, double*, int)`|
|124|0.00|100.05|0.00|1|`perf::perf()`|
|125|0.00|100.05|0.00|1|`perf::~perf()`|
|126|0.00|100.05|0.00|1|`Ansorg::set_ABp()`|
|127|0.00|100.05|0.00|1|`Ansorg::Ansorg(char*, int)`|
|128|0.00|100.05|0.00|1|`Ansorg::~Ansorg()`|
|129|0.00|100.05|0.00|1|`monitor::print_message(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)`|
|130|0.00|100.05|0.00|1|`std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, int, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, int> > >::~map()`|
|131|0.00|100.05|0.00|1|`std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >::pair<char const (&) [9], std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, true>(char const (&) [9], std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&)`|
