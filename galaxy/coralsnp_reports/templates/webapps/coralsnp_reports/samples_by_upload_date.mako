<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">Samples</h3>
        <h4 align="center">
            Click the number of samples to view
            the number of samples uploaded per month
        </h4>
        <table align="center" class="colored">
            %if num_samples == 0:
                <tr><td>There are no uploaded samples</td></tr>
            %else:
                <tr class="header">
                    <td align="center">
                        Number of Samples
                    </td>
                </tr>
                <tr class="tr">
                    <td align="center">
                        <a href="${h.url_for( controller='samples', action='per_month', sort_id='default', order='default' )}">
                            ${num_samples}
                        </a>
                    </td>
                </tr>
            %endif
        </table>
    </div>
</div>

