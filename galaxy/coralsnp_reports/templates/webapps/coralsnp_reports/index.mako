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
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='samples', action='all', sort_id='default', order='default' )}">All Samples</a></div>
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='samples', action='specified_date', specified_date=datetime.utcnow().strftime( "%Y-%m-%d" ), sort_id='default', order='default' )}">Today's Samples</a></div>
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='samples', action='specified_month', sort_id='default', order='default' )}">Samples per day this month</a></div>
                        <div class="toolTitle"><a target="galaxy_main" href="${h.url_for( controller='samples', action='per_month', sort_id='default', order='default' )}">Samples per month</a></div>
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
