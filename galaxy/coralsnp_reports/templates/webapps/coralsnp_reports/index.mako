<%inherit file="/webapps/coralsnp_reports/base_panels.mako"/>

<%def name="init()">
    ${parent.init()}
    <%
        self.has_left_panel=True
        self.has_right_panel=False
        self.active_view="coralsnp_reports"
    %>
</%def>

<%def name="stylesheets()">
    ${parent.stylesheets()}
    ## Include "base.css" for styling tool menu and forms (details)
    ${h.css( "base" )}

    ## But make sure styles for the layout take precedence
    ${parent.stylesheets()}

</%def>

<%def name="javascripts()">
    ${parent.javascripts()}
</%def>

<%def name="left_panel()">
    <%
        from datetime import datetime
        from time import mktime, strftime, localtime
    %>
    <div class="unified-panel-header" unselectable="on">
        <div class='unified-panel-header-inner'><span>Reports</span>
            <a target="galaxy_main" href="${h.url_for( controller='root', action='index' )}" class="float-right">
                <span class="fa fa-home"></span>
            </a>
        </div>
    </div>
    <div class="unified-panel-body">
        <div class="toolMenu">
            <div class="toolSectionList">
                <div class="toolSectionPad"></div>
                <div class="toolSectionTitle">
                    <span>Samples</span>
                </div>
                <div class="toolSectionBody">
                    <div class="toolSectionBg">
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='samples', action='all', sort_id='default', order='default' )}">All uploaded samples</a></div>
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='samples', action='all_by_upload_date', sort_id='default', order='default' )}">All samples by upload date</a></div>
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='samples', action='per_month', sort_id='default', order='default' )}">Samples uploaded per month</a></div>
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='samples', action='specified_month', sort_id='default', order='default' )}">Samples uploaded per day this month</a></div>
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='samples', action='specified_date', specified_date=datetime.utcnow().strftime( "%Y-%m-%d" ), sort_id='default', order='default' )}">Samples uploaded today</a></div>
                    </div>
                </div>
                <div class="toolSectionPad"></div>
                <div class="toolSectionTitle">
                    <span>Genotypes</span>
                </div>
                <div class="toolSectionBody">
                    <div class="toolSectionBg">
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='genotypes', action='all', sort_id='default', order='default' )}">All genotypes of uploaded samples</a></div>
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='genotypes', action='all_by_sample_upload_date', sort_id='default', order='default' )}">All genotypes by sample upload date</a></div>
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='genotypes', action='per_month', sort_id='default', order='default' )}">Genotypes by sample uploads per month</a></div>
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='genotypes', action='specified_month', sort_id='default', order='default' )}">Genotypes by sample uploads per day this month</a></div>
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='genotypes', action='specified_date', specified_date=datetime.utcnow().strftime( "%Y-%m-%d" ), sort_id='default', order='default' )}">Genotypes for samples uploaded today</a></div>
                    </div>
                </div>
                <div class="toolSectionPad"></div>
                <div class="toolSectionTitle">
                    <span>Phenotypes</span>
                </div>
                <div class="toolSectionBody">
                    <div class="toolSectionBg">
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='phenotypes', action='all', sort_id='default', order='default' )}">All phenotypes of uploaded samples</a></div>
                    </div>
                </div>
                <div class="toolSectionPad"></div>
                <div class="toolSectionTitle">
                    <span>Collectors</span>
                </div>
                <div class="toolSectionBody">
                    <div class="toolSectionBg">
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for(controller='collectors', action='all', sort_id='default', order='default')}">All collectors</a></div>
                    </div>
                </div>
                <div class="toolSectionPad"></div>
                <div class="toolSectionTitle">
                    <span>Experiments</span>
                </div>
                <div class="toolSectionBody">
                    <div class="toolSectionBg">
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='experiments', action='all', sort_id='default', order='default' )}">All experiments for uploaded samples</a></div>
                    </div>
                </div>
                <div class="toolSectionPad"></div>
                <div class="toolSectionTitle">
                    <span>Colonies</span>
                </div>
                <div class="toolSectionBody">
                    <div class="toolSectionBg">
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='colonies', action='all', sort_id='default', order='default' )}">All colonies of uploaded samples</a></div>
                    </div>
                </div>
                <div class="toolSectionPad"></div>
                <div class="toolSectionTitle">
                    <span>Reefs</span>
                </div>
                <div class="toolSectionBody">
                    <div class="toolSectionBg">
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='reefs', action='all', sort_id='default', order='default' )}">All reefs of uploaded samples</a></div>
                    </div>
                </div>
                <div class="toolSectionPad"></div>
                <div class="toolSectionTitle">
                    <span>Taxonomies</span>
                </div>
                <div class="toolSectionBody">
                    <div class="toolSectionBg">
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='taxonomies', action='all', sort_id='default', order='default' )}">All taxonomies of uploaded samples</a></div>
                    </div>
                </div>
            </div>
        </div>
    </div>
</%def>

<%def name="center_panel()">
    <% center_url = h.url_for( controller='root', action='home' ) %>
    <iframe name="galaxy_main" id="galaxy_main" frameborder="0" style="position: absolute; width: 100%; height: 100%;" src="${center_url}"> </iframe>
</%def>

