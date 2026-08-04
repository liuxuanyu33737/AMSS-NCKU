# AMSS-NCKU 8个MPI进程 gprof 综合报告

- 成功分析的MPI进程数：8
- 综合统计函数数量：131
- 运行配置：8个MPI进程、20个演化步
- 分析对象：CPU程序 `ABE`

## 综合函数耗时排名

|排名|总体占比(%)|累计占比(%)|平均占比(%)|最小占比(%)|最大占比(%)|自身时间合计(s)|函数名称|
|---:|---:|---:|---:|---:|---:|---:|---|
|1|47.47|47.47|49.09|38.76|57.61|1813.98|`compute_rhs_bssn_`|
|2|14.95|62.41|12.52|0.00|27.97|571.16|`polint_`|
|3|7.27|69.69|7.54|5.93|9.02|277.88|`lopsided_`|
|4|5.91|75.60|6.12|4.78|7.30|225.86|`fdderivs_`|
|5|4.93|80.53|5.11|3.97|6.15|188.55|`kodis_`|
|6|4.82|85.35|4.85|4.60|5.17|184.06|`prolong3_`|
|7|4.38|89.72|4.54|3.53|5.40|167.27|`fderivs_`|
|8|1.89|91.61|1.94|1.55|2.27|72.04|`enforce_ga_`|
|9|1.76|93.37|1.83|1.37|2.25|67.16|`rungekutta4_rout_`|
|10|1.54|94.91|1.58|1.26|1.80|58.85|`symmetry_bd_`|
|11|1.43|96.34|1.32|0.73|2.02|54.75|`_init`|
|12|0.75|97.09|0.77|0.62|0.87|28.78|`restrict3_`|
|13|0.72|97.81|0.61|0.00|1.34|27.61|`polin3_`|
|14|0.46|98.27|0.47|0.38|0.56|17.49|`average2_`|
|15|0.30|98.57|0.25|0.00|0.56|11.47|`decide3d_`|
|16|0.28|98.85|0.29|0.22|0.44|10.72|`misc::fact(int)`|
|17|0.21|99.06|0.21|0.18|0.26|8.01|`admmass_bssn_`|
|18|0.17|99.23|0.17|0.12|0.20|6.53|`copy_`|
|19|0.16|99.39|0.17|0.13|0.20|6.15|`getnp4_`|
|20|0.11|99.50|0.10|0.08|0.14|4.06|`Patch::Interp_Points(MyList<var>*, int, double**, double*, int)`|
|21|0.08|99.58|0.07|0.00|0.17|3.19|`global_interp_`|
|22|0.07|99.66|0.07|0.05|0.11|2.73|`misc::Wigner_d_function(int, int, int, double)`|
|23|0.07|99.72|0.07|0.04|0.11|2.50|`surface_integral::surf_Wave(double, int, cgh*, var*, var*, int, int, int, double*, double*, monitor*)`|
|24|0.06|99.78|0.06|0.05|0.07|2.32|`Ansorg::interpolate_tri_bar(double, double, double, int, int, int, double*, double*, double*, double*)`|
|25|0.05|99.83|0.05|0.04|0.08|2.01|`lowerboundset_`|
|26|0.02|99.86|0.03|0.01|0.04|0.90|`sommerfeld_rout_`|
|27|0.02|99.88|0.02|0.01|0.03|0.88|`monitor::monitor(char const*, int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)`|
|28|0.02|99.90|0.02|0.01|0.03|0.82|`Parallel::transfer(MyList<Parallel::gridseg>**, MyList<Parallel::gridseg>**, MyList<var>*, MyList<var>*, int)`|
|29|0.02|99.92|0.02|0.01|0.03|0.64|`bssn_class::Step(int, int)`|
|30|0.01|99.93|0.01|0.01|0.02|0.53|`l2normhelper_`|
|31|0.01|99.95|0.02|0.00|0.03|0.53|`sommerfeld_routbam_`|
|32|0.01|99.96|0.01|0.01|0.02|0.47|`surface_integral::surf_MassPAng(double, int, cgh*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, var*, double*, monitor*)`|
|33|0.01|99.97|0.01|0.01|0.01|0.32|`get_ansorg_nbhs_`|
|34|0.00|99.97|0.01|0.00|0.01|0.17|`average_`|
|35|0.00|99.98|0.01|0.00|0.01|0.14|`Parallel::build_gstl(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>**, MyList<Parallel::gridseg>**)`|
|36|0.00|99.98|0.00|0.00|0.01|0.12|`Block::getdX(int)`|
|37|0.00|99.98|0.00|0.00|0.01|0.10|`Parallel::Sync(Patch*, MyList<var>*, int)`|
|38|0.00|99.98|0.00|0.00|0.01|0.09|`Parallel::prepare_inter_time_level(Patch*, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, int)`|
|39|0.00|99.99|0.00|0.00|0.01|0.08|`Parallel::gs_subtract(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*)`|
|40|0.00|99.99|0.00|0.00|0.01|0.07|`Ansorg::find_point_bisection(double, int, double*, int)`|
|41|0.00|99.99|0.00|0.00|0.01|0.07|`Patch::Interp_ONE_Point(MyList<var>*, double*, double*, int)`|
|42|0.00|99.99|0.00|0.00|0.00|0.04|`Block::swapList(MyList<var>*, MyList<var>*, int)`|
|43|0.00|99.99|0.00|0.00|0.01|0.04|`Parallel::build_ghost_gsl(Patch*)`|
|44|0.00|99.99|0.00|0.00|0.01|0.04|`Parallel::build_owned_gsl2(Patch*, int)`|
|45|0.00|99.99|0.00|0.00|0.00|0.03|`Parallel::Restrict(MyList<Patch>*, MyList<Patch>*, MyList<var>*, MyList<var>*, int)`|
|46|0.00|99.99|0.00|0.00|0.00|0.03|`Parallel::Sync(MyList<Patch>*, MyList<var>*, int)`|
|47|0.00|100.00|0.00|0.00|0.00|0.03|`cgh::Interp_One_Point(MyList<var>*, double*, double*, int)`|
|48|0.00|100.00|0.00|0.00|0.00|0.02|`Ansorg::ps_u_at_xyz(double, double, double)`|
|49|0.00|100.00|0.00|0.00|0.00|0.02|`Parallel::build_bulk_gsl(Block*, Patch*)`|
|50|0.00|100.00|0.00|0.00|0.00|0.02|`Parallel::build_complete_gsl(Patch*)`|
|51|0.00|100.00|0.00|0.00|0.01|0.02|`Parallel::gsl_subtract(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*)`|
|52|0.00|100.00|0.00|0.00|0.00|0.02|`bssn_class::Constraint_Out()`|
|53|0.00|100.00|0.00|0.00|0.00|0.01|`Ansorg::xyz_to_ABp(double, double, double, double*, double*, double*)`|
|54|0.00|100.00|0.00|0.00|0.00|0.01|`Parallel::Dump_Data(Patch*, MyList<var>*, char*, double, double, int)`|
|55|0.00|100.00|0.00|0.00|0.00|0.01|`Parallel::build_buffer_gsl(Patch*)`|
|56|0.00|100.00|0.00|0.00|0.00|0.01|`Parallel::build_owned_gsl0(Patch*, int)`|
|57|0.00|100.00|0.00|0.00|0.00|0.01|`Parallel::build_owned_gsl4(Patch*, int, int)`|
|58|0.00|100.00|0.00|0.00|0.00|0.01|`Patch::getdX(int)`|
|59|0.00|100.00|0.00|0.00|0.00|0.01|`bssn_class::Compute_Constraint()`|
|60|0.00|100.00|0.00|0.00|0.00|0.01|`bssn_class::RecursiveStep(int)`|
|61|0.00|100.00|0.00|0.00|0.00|0.00|`Ansorg::Ansorg(char*, int)`|
|62|0.00|100.00|0.00|0.00|0.00|0.00|`Ansorg::set_ABp()`|
|63|0.00|100.00|0.00|0.00|0.00|0.00|`Ansorg::~Ansorg()`|
|64|0.00|100.00|0.00|0.00|0.00|0.00|`Block::Block(int, int*, double*, int, int, int, int, int)`|
|65|0.00|100.00|0.00|0.00|0.00|0.00|`Block::~Block()`|
|66|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::Dump_Data(MyList<Patch>*, MyList<var>*, char*, double, double)`|
|67|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::KillBlocks(MyList<Patch>*)`|
|68|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::L2Norm(Patch*, var*)`|
|69|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::OutBdLow2Hi(Patch*, Patch*, MyList<var>*, MyList<var>*, int)`|
|70|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::PatList_Interp_Points(MyList<Patch>*, MyList<var>*, int, double**, double*, int)`|
|71|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::add_ghost_touch(MyList<Parallel::gridseg>*&)`|
|72|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::aligncheck(double*, double*, int, double*, int*)`|
|73|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::build_complete_gsl_virtual(MyList<Patch>*)`|
|74|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::checkgsl(MyList<Parallel::gridseg>*, bool)`|
|75|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::cut_gs(MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*, MyList<Parallel::gridseg>*&)`|
|76|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::cut_gsl(MyList<Parallel::gridseg>*&)`|
|77|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::distribute(MyList<Patch>*, int, int, int, bool, int)`|
|78|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::fill_level_data(MyList<Patch>*, MyList<Patch>*, MyList<Patch>*, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, int, bool, bool)`|
|79|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::merge_gsl(MyList<Parallel::gridseg>*&, double)`|
|80|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::partition3(int*, int, int*, int, int*)`|
|81|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::prepare_inter_time_level(Patch*, MyList<var>*, MyList<var>*, MyList<var>*, int)`|
|82|0.00|100.00|0.00|0.00|0.00|0.00|`Parallel::writefile(double, int, int, int, double, double, double, double, double, double, char*, double*)`|
|83|0.00|100.00|0.00|0.00|0.00|0.00|`Patch::Patch(int, int*, double*, int, bool, int)`|
|84|0.00|100.00|0.00|0.00|0.00|0.00|`Patch::checkPatch(bool)`|
|85|0.00|100.00|0.00|0.00|0.00|0.00|`Patch::~Patch()`|
|86|0.00|100.00|0.00|0.00|0.00|0.00|`bssn_class::Compute_Psi4(int)`|
|87|0.00|100.00|0.00|0.00|0.00|0.00|`bssn_class::Evolve(int)`|
|88|0.00|100.00|0.00|0.00|0.00|0.00|`bssn_class::Initialize()`|
|89|0.00|100.00|0.00|0.00|0.00|0.00|`bssn_class::Interp_Constraint(bool)`|
|90|0.00|100.00|0.00|0.00|0.00|0.00|`bssn_class::Read_Ansorg()`|
|91|0.00|100.00|0.00|0.00|0.00|0.00|`bssn_class::Setup_Black_Hole_position()`|
|92|0.00|100.00|0.00|0.00|0.00|0.00|`bssn_class::bssn_class(double, double, double, double, double, double, double, int, int, char*, double, double, double, int, int, int, double, double)`|
|93|0.00|100.00|0.00|0.00|0.00|0.00|`bssn_class::check_Stdin_Abort()`|
|94|0.00|100.00|0.00|0.00|0.00|0.00|`bssn_class::compute_Porg_rhs(double**, double**, var*, var*, var*, int)`|
|95|0.00|100.00|0.00|0.00|0.00|0.00|`bssn_class::~bssn_class()`|
|96|0.00|100.00|0.00|0.00|0.00|0.00|`cgh::Regrid_Onelevel(int, int, int, double**, double**, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, bool, monitor*)`|
|97|0.00|100.00|0.00|0.00|0.00|0.00|`cgh::cgh(int, int, int, char*, int, monitor*)`|
|98|0.00|100.00|0.00|0.00|0.00|0.00|`cgh::compose_cgh(int)`|
|99|0.00|100.00|0.00|0.00|0.00|0.00|`cgh::construct_patchlist(int, int)`|
|100|0.00|100.00|0.00|0.00|0.00|0.00|`cgh::read_bbox(int, char*)`|
|101|0.00|100.00|0.00|0.00|0.00|0.00|`cgh::recompose_cgh_Onelevel(int, int, MyList<var>*, MyList<var>*, MyList<var>*, MyList<var>*, int, bool)`|
|102|0.00|100.00|0.00|0.00|0.00|0.00|`cgh::sethandle(monitor*)`|
|103|0.00|100.00|0.00|0.00|0.00|0.00|`cgh::settrfls(int)`|
|104|0.00|100.00|0.00|0.00|0.00|0.00|`cgh::~cgh()`|
|105|0.00|100.00|0.00|0.00|0.00|0.00|`checkpoint::addvariablelist(MyList<var>*)`|
|106|0.00|100.00|0.00|0.00|0.00|0.00|`checkpoint::checkpoint(bool, char const*, int)`|
|107|0.00|100.00|0.00|0.00|0.00|0.00|`checkpoint::~checkpoint()`|
|108|0.00|100.00|0.00|0.00|0.00|0.00|`frame_dummy`|
|109|0.00|100.00|0.00|0.00|0.00|0.00|`misc::gaulegf(double, double, double*, double*, int)`|
|110|0.00|100.00|0.00|0.00|0.00|0.00|`misc::inversearray(double*, int)`|
|111|0.00|100.00|0.00|0.00|0.00|0.00|`misc::parse_parts(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, int&)`|
|112|0.00|100.00|0.00|0.00|0.00|0.00|`misc::parse_parts(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, int&, int&, int&)`|
|113|0.00|100.00|0.00|0.00|0.00|0.00|`monitor::print_message(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)`|
|114|0.00|100.00|0.00|0.00|0.00|0.00|`monitor::writefile(double, int, double*)`|
|115|0.00|100.00|0.00|0.00|0.00|0.00|`monitor::writefile(double, int, double*, double*)`|
|116|0.00|100.00|0.00|0.00|0.00|0.00|`monitor::~monitor()`|
|117|0.00|100.00|0.00|0.00|0.00|0.00|`perf::MemoryUsage(unsigned long*, unsigned long*, unsigned long*, unsigned long*, unsigned long*, unsigned long*, int)`|
|118|0.00|100.00|0.00|0.00|0.00|0.00|`perf::perf()`|
|119|0.00|100.00|0.00|0.00|0.00|0.00|`perf::sample_mem_usage(int)`|
|120|0.00|100.00|0.00|0.00|0.00|0.00|`perf::~perf()`|
|121|0.00|100.00|0.00|0.00|0.00|0.00|`setpbh(int, double**, double*, int)`|
|122|0.00|100.00|0.00|0.00|0.00|0.00|`std::_Rb_tree<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::_Select1st<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >::find(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&)`|
|123|0.00|100.00|0.00|0.00|0.00|0.00|`std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, int, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, int> > >::~map()`|
|124|0.00|100.00|0.00|0.00|0.00|0.00|`std::pair<std::_Rb_tree_iterator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, bool> std::_Rb_tree<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::_Select1st<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >::_M_insert_unique<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >(std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >&&)`|
|125|0.00|100.00|0.00|0.00|0.00|0.00|`std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >::pair<char const (&) [9], std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&, true>(char const (&) [9], std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >&)`|
|126|0.00|100.00|0.00|0.00|0.00|0.00|`std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >::~pair()`|
|127|0.00|100.00|0.00|0.00|0.00|0.00|`surface_integral::surface_integral(int)`|
|128|0.00|100.00|0.00|0.00|0.00|0.00|`surface_integral::~surface_integral()`|
|129|0.00|100.00|0.00|0.00|0.00|0.00|`var::setpropspeed(double)`|
|130|0.00|100.00|0.00|0.00|0.00|0.00|`var::var(char const*, int, double, double, double)`|
|131|0.00|100.00|0.00|0.00|0.00|0.00|`var::~var()`|

## 指标说明

- **总体占比**：该函数在8个进程中的自身时间合计，占全部8个进程采样时间的比例。
- **平均占比**：该函数在8个MPI进程中占比的算术平均值。
- **最小/最大占比**：用于观察不同MPI进程之间是否存在负载差异。
